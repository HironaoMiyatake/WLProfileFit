import corner
import numpy as np
import sys
import matplotlib.pyplot as plt
import os

# input
filename = sys.argv[1]

# load data
param_names = [r"$M_{200m}$", r"$c_{200m}$"]
chain = np.load(filename)
samples = chain.reshape((-1, chain.shape[2]))

# plot triangle
fig = corner.corner(samples, labels = param_names)

# calculate statistics
samples[:, 0] = samples[:, 0]
M_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
print "M [10^14 Msun/h]:", M_mcmc
print "c:", c_mcmc

filename = os.path.join(os.path.dirname(filename), "contour.png")
plt.savefig(filename)
print "file written in:", filename

#plt.show()
