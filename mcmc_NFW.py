import os, sys
import numpy as np
import scipy.interpolate as spInterpolate
import NFW
import emcee
from emcee.utils import MPIPool

# input file
filename = sys.argv[1]

# configurations
rmin = 0.2
rmax = 10.0
# directory name prefix
dir_prefix = "fit_NFW"

# chain will be seved every nstep. In total nbunch * nstep samplings.
nbunch = 10
nstep = 50
nwalkers = 96

# initial value for parameters
init_M = 2.0 #[10^14 Msun/h]
init_c = 4.3

# main body
# definition of model
nfw = NFW.NFW(200, 0.3, 0.7)
def getModel(x, M, c):
    return nfw.deltaSigma(x, M*10**14, c)

# likelihood
def lnlike(theta, x, y, yerr):
    M, c = theta
    model = getModel(x, M, c)
    diff = y-model
    return -0.5*np.sum(diff**2/yerr**2)

# prior
def lnprior(theta):
    M, c = theta
    if 0.05 < M < 10 and 0.5 < c < 7.0:
        return 0.0
    return -np.inf

# posterior
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    lpb = lp + lnlike(theta, x, y, yerr)
    return lpb

pool = MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)

# read data
data = np.genfromtxt(filename, names = True)
sel = (data["r"]> rmin) & (data["r"] < rmax)
x = data["r"][sel]
y = data["dSigma"][sel]
yerr = data["dSigma_err"][sel]

# setup output directory
dir = "%s_%s" % (dir_prefix, filename)
print "output directory: %s" % dir
if os.path.exists(dir) == False:
    os.mkdir(dir)

# write down data
np.savetxt(os.path.join(dir, "data.dat"), np.array([x,y, yerr]).transpose())

# run MCMC
ndim = 2
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr), pool = pool)
for i in range(nbunch):
    if i == 0:
        pos = [[init_M, init_c] + 1e-4*np.random.randn(ndim) for j in range(nwalkers)]
    else :
        pos = sampler.chain[:,-1,:]
    sampler.run_mcmc(pos, nstep)
    filename_bunch_chains = os.path.join(dir, "chains")
    np.save(filename_bunch_chains, sampler.chain)
    filename_bunch_lnprobabilities = os.path.join(dir, "lnprobabilities")
    np.save(filename_bunch_lnprobabilities, sampler.lnprobability)
    print "%s/%s bunch completed. File written in %s.npy" % (i+1, nbunch, filename_bunch_chains)
pool.close()
