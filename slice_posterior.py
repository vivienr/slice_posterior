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

parser.add_argument('-q', type=float, help='mass ratio m_2/m_1 < 1 (SMALLER THAN 1).')
parser.add_argument('-a1z', type=float, help='a1z of first BH (BETWEEN -1 AND 1)')
parser.add_argument('-a2z', type=float, help='a2z of second BH (BETWEEN -1 AND 1)')
parser.add_argument('--tilt1', type=float, default=None, help='BH 1 tilt-angle used by the NR simulation.  If not specified, then this dimension is not sliced')
parser.add_argument('--tilt2', type=float, default=None, help='BH 2 tilt-angle used by the NR simulation.  If not specified, then this dimension is not sliced')
parser.add_argument('--tilt_factor', type=float, default=2.,
                    help="allow range of tilt1/2 to be this factor larger than the range in 1/q, a1z and a2z.  Default=2.")

parser.add_argument('-c', '--credible-interval', type=float, help='Credible interval to use as subset of the posterior.',
                    default=0.9)
parser.add_argument('-d', '--delta', type=float, default=None,
                    help='Allowed +/- range around the given mass-ratio, spin1, spin2 values.')
parser.add_argument('--Nremaining', type=int, default=None,
                    help='Choose a delta such that Nremaining posterior samples remain in sliced posterior')

args = parser.parse_args()

# parse Nremaining and delta
if args.delta is not None:
    delta=args.d
    if args.Nremaining is not None:
        raise Exception('Must not specify both --delta and --Nremaining')
    Nremaining=None
else:
    if args.Nremaining is None:
        raise Exception('Must specify one of {--delta,  --Nremaining}')
    Nremaining=args.Nremaining
    delta=None

# shortcuts for other parameters
q=args.q
if(q>1.):
    raise Exception("Mass-ratio -q must be <=1.  q=m2/m1, with m1 being the more massive BH")
a1z=args.a1z
a2z=args.a2z
tilt1=args.tilt1
tilt2=args.tilt2
tilt_factor=args.tilt_factor


################################################################
# function to compute sliced posterior
def Do_Slicing(q, a1z, a2z, tilt1, tilt2, delta, tilt_factor, data):
    condition=[np.abs(q/data['q']-1.)<delta,np.abs(data['a1z']-a1z)<delta,np.abs(data['a2z']-a2z)<delta]

    if(tilt1 is not None):
        condition.append(np.abs(data['tilt1']-tilt1)<delta*tilt_factor)
    if(tilt2 is not None):
        condition.append(np.abs(data['tilt2']-tilt2)<delta*tilt_factor)
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
    out=Do_Slicing(q, a1z, a2z, tilt1, tilt2, delta, tilt_factor, data)
else:
    print "Determine delta by bisection.  Goal Nremaining={}".format(Nremaining)
    # need to bisect
    delta_min=0.
    delta_max=1.
    out_min=Do_Slicing(q, a1z, a2z, tilt1, tilt2, delta_min, tilt_factor, data)
    out_max=Do_Slicing(q, a1z, a2z, tilt1, tilt2, delta_max, tilt_factor, data)
    if len(out_min) > Nremaining:
        raise Exception('logic error in algorithm')
    if len(out_max) < Nremaining:
        raise Exception('logic error in algorithm')
    while(len(out_max) != len(out_min)):
        print "delta-range=[{:.6f}, {:.6f}], Nsamples=[{}, {}]".format(delta_min, delta_max, len(out_min), len(out_max))

        # bifurcate
        delta=(delta_min+delta_max)/2.
        out=Do_Slicing(q,a1z, a2z, tilt1, tilt2, delta, tilt_factor, data)
        if len(out)==Nremaining: break
        if len(out)<Nremaining:
           delta_min=delta
           out_min=out
        else:
           delta_max=delta
           out_max=out


if len(out) == 0:
    print "No sample point found within +/- "+str(delta)+" of q="+str(1./q)+", a1z="+str(a1z)+", a2z="+str(a2z)
    os.sys.exit()
else:
    print "Highest posterior point found within +/- "+str(delta)+" of q="+str(q)+", a1z="+str(a1z)+", a2z="+str(a2z)+":"
    if(tilt1 is not None) or (tilt2 is not None):
        print "['q', 'a1z', 'a2z', 'tilt1', 'tilt2']={}".format(out[['q', 'a1z', 'a2z', 'tilt1', 'tilt2']][0])
    else:
        print "['q', 'a1z', 'a2z']={}".format(out[['q', 'a1z', 'a2z']][0])
    print "['mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi']={}".format(out[['mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi']][0])
    print "['time','h1_end_time','l1_end_time','h1l1_delay']={}".format(out[['time','h1_end_time','l1_end_time','h1l1_delay']][0])

    maxLpost=(data['logprior']+data['logl'])[0]
    maxLslice=(out['logprior']+out['logl'])[0]
    print "overall highest logposterior ={}".format(maxLpost)
    print "highest logposterior in slice={}".format(maxLslice)
    print "difference in logposterior={}".format(maxLpost-maxLslice)
    
    #    print out[['q', 'a1z', 'a2z', 'mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi', 'time','h1_end_time','l1_end_time','h1l1_delay']][0]

if args.output:
    print "Writing "+str(len(out))+" points to file "+args.output
    np.savetxt(args.output,out,header='\t'.join(out.dtype.names))

pts=np.array([[x,y,z] for x,y,z in cred_data[['q','a1z','a2z']]])
mtot=griddata(pts, cred_data['mtotal'], (q, a1z, a2z), method='linear')

print "Interpolated total mass: "+str(mtot)
