#!/usr/bin/env python

from __future__ import division

import numpy as np
from scipy.interpolate import griddata
import argparse
import os

parser = argparse.ArgumentParser(description=
                                 """Slices and interpolate a LIGO CBC-PE posterior file to give the
parameters (total mass and extrinsic parameters) fitting best specific
mass-ratio, a1z, a2z, and (if desired) tilt1 and tilt2 values.
""")



parser.add_argument('-p', '--posterior', type=str, help='posterior file (ASCII table with one-line header).')

parser.add_argument('-o', '--output', type=str, nargs='?', help='output file.',
                    default=None)

parser.add_argument('-c', '--credible-interval', type=float, help='Credible interval to use as subset of the posterior.',
                    default=0.9)
group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('--delta', type=float, default=None,
                   help='Allowed +/- range around the given mass-ratio, spin1, spin2 values.')
group.add_argument('--Nremaining', type=int, default=None,
                    help='Choose a delta such that Nremaining posterior samples remain in sliced posterior')

# list of slice-able variables.  This list
# can be extended as desired

VARS=['q', 'a1z', 'a2z', 'tilt1', 'tilt2']
for v in VARS:
    parser.add_argument('-{}'.format(v), type=float, default=None,
                        help="specify a value for {} at which to slice".format(v))
for v in VARS:
    parser.add_argument('-frac_{}'.format(v), default=False, action='store_true',
                        help='enforce a fractional width for variable {}'.format(v))
for v in VARS:
    parser.add_argument('-mult_{}'.format(v), default=1., type=float,
                        help="width multiplier of this variable in the posterior slice (default: %(default)s)")


args = parser.parse_args()

delta=args.delta
Nremaining=args.Nremaining
argsd=vars(args) # turn args into dictionary for use in Do_Slicing


################################################################
# function to compute sliced posterior
def Do_Slicing(data, delta):
    condition=[]
    for v in VARS:
        if argsd[v] is not None:
            if argsd['frac_'+v]:
                # fractional width
                condition.append(np.abs((argsd[v]-data[v])/argsd[v])<delta*argsd['mult_'+v])
            else:
                # absolute width
                condition.append(np.abs(argsd[v]-data[v])<delta*argsd['mult_'+v])
    out=data[np.where(np.all(condition,axis=0))]
    return out
################################################################


################################################################
# Import data
org_data=np.genfromtxt(args.posterior, names=True)

# Sort the data in order of increasing logposterior=logprior+loglikelihood
sorted_data=org_data[np.argsort(org_data['logprior']+org_data['logl'])][::-1]
# Take the top (credible-interval)
cred_data=sorted_data[0:int(args.credible_interval*len(sorted_data))]

data=cred_data

################################################################
# do slicing
if delta is not None:
    out=Do_Slicing(data, delta)
else:
    print "Determine delta by bisection.  Goal Nremaining={}".format(Nremaining)
    # need to bisect
    delta_min=0.
    delta_max=1.
    out_min=Do_Slicing(data, delta_min)
    out_max=Do_Slicing(data, delta_max)
    if len(out_min) > Nremaining:
        raise Exception('logic error in algorithm')
    if len(out_max) < Nremaining:
        raise Exception('logic error in algorithm')
    while(len(out_max) != len(out_min)):
        print "delta-range=[{:.6f}, {:.6f}], Nsamples=[{}, {}]".format(delta_min, delta_max, len(out_min), len(out_max))

        # bifurcate
        delta=(delta_min+delta_max)/2.
        out=Do_Slicing(data, delta)
        if len(out)==Nremaining: break
        if len(out)<Nremaining:
           delta_min=delta
           out_min=out
        else:
           delta_max=delta
           out_max=out

################
# generate nice output
s=""
for v in VARS:
    if argsd[v] is None: continue
    range=delta*argsd['mult_'+v]
    if argsd['frac_'+v]: range=range*argsd[v]
    s=s+", {}={}+-{:.4f}".format(v, argsd[v],range)
s=s[2:]

def PrintOnlyExisting(data, members):
    q=[]
    for m in members:
        if m in data.dtype.names: q.append(m)
    print "{}={}".format(q,data[q][0])


if len(out) == 0:
    print "No sample point found within "+s
    os.sys.exit()
else:
    print "Highest posterior point found within "+s
    PrintOnlyExisting(data, ['q', 'a1z', 'a2z', 'tilt1', 'tilt2'])
    print "['mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi']={}".format(out[['mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi']][0])
    print "['time','h1_end_time','l1_end_time','h1l1_delay']={}".format(out[['time','h1_end_time','l1_end_time','h1l1_delay']][0])

    maxLpost=(data['logprior']+data['logl'])[0]
    maxLslice=(out['logprior']+out['logl'])[0]
    print "overall highest logposterior ={}".format(maxLpost)
    print "highest logposterior in slice={}".format(maxLslice)
    print "difference in logposterior={}".format(maxLpost-maxLslice)


if args.output:
    print "Writing "+str(len(out))+" points to file "+args.output
    np.savetxt(args.output,out,header='\t'.join(out.dtype.names))

#pts=np.array([[x,y,z] for x,y,z in cred_data[['q','a1z','a2z']]])
#mtot=griddata(pts, cred_data['mtotal'], (q, a1z, a2z), method='linear')
#print "Interpolated total mass: "+str(mtot)
