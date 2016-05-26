#!/usr/bin/env python

from __future__ import division

import numpy as np
from scipy.interpolate import griddata
import argparse

parser = argparse.ArgumentParser(description="Slices and interpolate a LIGO CBC-PE posterior file to give the parameters (total mass and extrinsic parameters) fitting best specific mass-ratio, spin1, spin2 values.")

parser.add_argument('-p', '--posterior', type=str, help='posterior file (ASCII table with one-line header).')

parser.add_argument('-q', type=float, help='mass ratio (GREATER THAN 1) used by the NR simulation.')
parser.add_argument('-chi1', type=float, help='spin magnitude (BETWEEN -1 AND 1) used by the NR simulation.')
parser.add_argument('-chi2', type=float, help='spin magnitude (BETWEEN -1 AND 1) used by the NR simulation.')

parser.add_argument('-c', '--credible-interval', type=float, help='Credible interval to use as subset of the posterior.',
                    default=0.9)
parser.add_argument('-d', '--delta', type=float, help='Allowed +/- range around the given mass-ratio, spin1, spin2 values.',
                    default=0.01)

parser.add_argument('-o', '--output', type=str, nargs='?', help='output file.',
                    default=None)

args = parser.parse_args()

org_data=np.genfromtxt(args.posterior, names=True)

# Sort the data in order of increasing logposterior=logprior+loglikelihood
sorted_data=org_data[np.argsort(org_data['logprior']+org_data['logl'])][::-1]
# Take the top (credible-interval)
cred_data=sorted_data[0:int(args.credible_interval*len(sorted_data))]

data=cred_data

delta=args.delta

q=1./args.q
chi1=args.chi1
chi2=args.chi2

condition=[np.abs(data['q']-q)<delta,np.abs(data['a1z']-chi1)<delta,np.abs(data['a2z']-chi2)<delta]

out=data[np.where(np.all(condition,axis=0))]

if len(out) == 0:
    print "No sample point found within +/- "+str(delta)+" of q="+str(1./q)+", chi1="+str(chi1)+", chi2="+str(chi2)
    exit()
else:
    print "Highest posterior point found within +/- "+str(delta)+" of q="+str(1./q)+", chi1="+str(chi1)+", chi2="+str(chi2)+":"
    print ['q', 'a1z', 'a2z', 'mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi', 'phi_orb', 'time','h1_end_time','l1_end_time','h1l1_delay']
    print out[['q', 'a1z', 'a2z', 'mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi', 'phi_orb', 'time','h1_end_time','l1_end_time','h1l1_delay']][0]

if args.output:
    print "Writing "+str(len(out))+" points to file "+args.output
    np.savetxt(args.output,out,header='\t'.join(out.dtype.names))

pts=np.array([[x,y,z] for x,y,z in cred_data[['q','a1z','a2z']]])
mtot=griddata(pts, cred_data['mtotal'], (q, chi1, chi2), method='linear')

print "Interpolated total mass: "+str(mtot)
