This is a code to fit a weak lesing profile. Currently only the NFW profile is supported.

You need scipy, emcee, and mpiexec to run this code.

# How to run this code?
You can learn by running the following test.
```bash
mpiexec -n 2 python mcmc_NFW.py test.dat
```
Chains are saved in `fit_NFW_test.dat`. Then show chains
```bash
python show_chains_NFW.py fit_NFW_test.dat/chains.npy
```
. You can cut out chains before burn-in by
```bash
python cutout_burn_in.py fit_NFW_test.dat 100
```
.