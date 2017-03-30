This is a code to fit a weak lesing profile. Currently only the NFW profile is supported.

You need scipy, emcee, and mpiexec to run this code.

# How to run this code?
You can learn by running the following test.
```bash
mpiexec -n 2 python mcmc_NFW.py test.dat
```
The format of input file should be
`r [Mpc/h] deltaSigma [h Msun/pc^2] deltaSigma_err [h Msun/pc^2]`
Chains are saved in `fit_NFW_test.dat`. Then show chains
```bash
python show_chains_NFW.py fit_NFW_test.dat/chains.npy
```
. You can cut out chains before burn-in by
```bash
python cutout_burn_in.py fit_NFW_test.dat 100
```
. In this case, the first 100 chains will be removed, and the output file is created at `fit_NFW_test.dat/chains.burnin100.npy`. You can then plot the contour by
```bash
python plot_triangle_NFW.py fit_NFW_test.dat/chains.burnin100.npy 
```
, and the fitting curve with the data points by
```bash
python plot_fit_NFW.py fit_NFW_test.dat/chains.burnin100.npy 
```
