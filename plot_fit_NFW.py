import os, sys
import numpy as np
import scipy.interpolate as spInterpolate
import NFW
import offcenteredNFW
import matplotlib.pyplot as plt
import re

# input
#filename_data = sys.argv[1]
filename_mcmc = sys.argv[1]
filename_data = os.path.join(os.path.dirname(filename_mcmc), "data.dat")
filename_lnp = re.sub("chains", "lnprobabilities", filename_mcmc)

# load data
data = np.genfromtxt(filename_data, names = ("r", "dSigma", "dSigma_err"))
x = data["r"]
y = data["dSigma"]
yerr = data["dSigma_err"]

# load chains
param_names = [r"$M_{200m}$", r"$c_{200m}$"]
chain = np.load(filename_mcmc)
samples = chain.reshape((-1, chain.shape[2]))
M = samples[:,0]
c = samples[:,1]

lnprobability = np.load(filename_lnp)
lnprobabilities = lnprobability.flatten()
lnp_argmax = np.argmax(lnprobabilities)
chi2 = -2.*lnprobabilities[lnp_argmax]
print "minimum chi2 from chains:", chi2

M_med = np.percentile(M, 50.)
c_med = np.percentile(c, 50.)
M_err, c_err = map(lambda v: (v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))

M_best = samples[lnp_argmax,0]
c_best = samples[lnp_argmax,1]

# definition of model
nfw = NFW.NFW(200, 0.315, 0.685)

x_model = np.logspace(-2, 2, 100)
NFW_cen = nfw.deltaSigma(x_model, M_best*10**14, c_best)

plt.errorbar(x, y, yerr, fmt = "o")
plt.plot(x_model, NFW_cen, "g-")
plt.loglog(nonposy='clip')
plt.xlim(0.05, 50)
plt.ylim(0.1, 300)
plt.xlabel(r"$r\ [h^{-1}\mathrm{comoving\ Mpc}]$")
plt.ylabel(r"$\Delta \Sigma_+\ [h\ M_\odot/\mathrm{comoving\ pc}^2]$")

model_tot = nfw.deltaSigma(x, M_best*10**14, c_best)
diff = model_tot - y
#chi2 = np.dot(diff, np.dot(invcov, np.transpose(diff)))
chi2 = np.sum((model_tot - y)**2/yerr**2)
print "chi2/dof: %.2f/%d = %s" % (chi2, len(x)-2, chi2/(len(x)-2.))

plt.annotate(r"$\chi^2/dof = %.2f/%d$" % (chi2, len(x)-2), xy = (0.65, 0.95), xycoords = "axes fraction")
plt.annotate(r"$M_{200m} = %.2f^{+%.2f}_{-%.2f} 10^{14}M_\odot/h$" % (M_med, M_err[0], M_err[1]), xy = (0.65, 0.9), xycoords = "axes fraction")
plt.annotate(r"$c_{200m} = %.2f^{+%.2f}_{-%.2f}$" % (c_med, c_err[0], c_err[1]), xy = (0.65, 0.85), xycoords = "axes fraction")

filename = os.path.join(os.path.dirname(filename_mcmc), "data_model.png")
plt.savefig(filename)
print "file written in:", filename

filename = os.path.join(os.path.dirname(filename_mcmc), "model.dat")
np.savetxt(filename, np.array([x_model, NFW_cen]).transpose())
print "file written in:", filename

plt.show()




