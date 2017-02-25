# slice_posterior

## Example

compute the subset of input-posterior samples, which have a
some parameters limited to certain user-specified values.

The user may either specify the width of the slice in each dimension
using -delta and -mult_VAR options, or the user can specify the
desired number of samples within the slice, in which case the
script determines the appropriate delta.

Running:

`slice_posterior.py -q 0.7 -frac_q -mult_q=2 -a1z 0.5 -a2z -0.5 --Nremaining 100 -p aligned_prior.dat`

Will generate the output:

```
Determine delta by bisection.  Goal Nremaining=100
delta-range=[0.000000, 1.000000], Nsamples=[0, 3053]
delta-range=[0.000000, 0.500000], Nsamples=[0, 1349]
delta-range=[0.000000, 0.250000], Nsamples=[0, 251]
delta-range=[0.125000, 0.250000], Nsamples=[42, 251]
delta-range=[0.125000, 0.187500], Nsamples=[42, 115]
delta-range=[0.156250, 0.187500], Nsamples=[73, 115]
delta-range=[0.171875, 0.187500], Nsamples=[93, 115]
delta-range=[0.171875, 0.179688], Nsamples=[93, 101]
delta-range=[0.175781, 0.179688], Nsamples=[96, 101]
delta-range=[0.177734, 0.179688], Nsamples=[98, 101]
delta-range=[0.178711, 0.179688], Nsamples=[99, 101]
delta-range=[0.178711, 0.179199], Nsamples=[99, 101]
delta-range=[0.178955, 0.179199], Nsamples=[99, 101]
delta-range=[0.179077, 0.179199], Nsamples=[99, 101]
Highest posterior point found within q=0.7+-0.2508, a1z=0.5+-0.1791, a2z=-0.5+-0.1791
['q', 'a1z', 'a2z']=(0.034280682065908, -0.497623181140534, -0.368766012712112)
['mtotal', 'theta_jn', 'distance', 'ra', 'dec', 'psi']=(45.05265626690489, 2.2521250523738314, 1970.2219978701826, 1.226130535759332, 0.283763239055403, 1.350051123907946)
['time','h1_end_time','l1_end_time','h1l1_delay']=(0.002765335016427, -0.0027795484131621704, -0.010016148170761578, -0.007236599757599407)
overall highest logposterior =32.059205
highest logposterior in slice=30.447156
difference in logposterior=1.612049
```

And create the `output.dat` file, which can then be plotted using https://github.com/vivienr/plot_posterior

`plot_posterior/plot_posterior.py -p output.dat aligned_prior.dat -n q a1z a2z -o plot.png`
![plot.png.png](https://github.com/vivienr/slice_posterior/blob/master/example/plot.png)

(see the example folder for `aligned_prior.dat` example input file)
