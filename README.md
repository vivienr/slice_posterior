# slice_posterior

## Example

Running:

`slice_posterior.py -p aligned_prior.dat -q 2 -chi1 0.5 -chi2 -0.5 -d 0.3 -c 1.0 -o output.dat`

Will generate the output:

```
Highest posterior point found within +/- 0.3 of q=2.0, chi1=0.5, chi2=-0.5:
['q', 'a1z', 'a2z', 'mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi', 'phi_orb', 'time', 'h1_end_time', 'l1_end_time', 'h1l1_delay']
(0.215562428819071, 0.209985184074858, -0.31854866419525, 35.58292842107642, 0.4241651829035867, 1983.4341437760845, 3.334889374419752, -0.5844587333721, 0.433760676566853, 3.726898854797899, -0.014173566088116, 0.005323766748355892, 0.007032821530383673, 0.0017090547820277809)
Writing 353 points to file output.dat
Interpolated total mass: 27.2920774447
```

And create the `output.dat` file

(see the example folder for `aligned_prior.dat` and `output.dat` data files)
