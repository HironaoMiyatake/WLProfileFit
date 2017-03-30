This is a code to fit a weak lesing profile. Currently only the NFW profile is supported.

You need scipy, emcee, and mpiexec to run this code.

# How to run this code?
You can learn by running the following test.
```bash
mpiexec -n 2 python mcmc_NFW.py test.dat
```