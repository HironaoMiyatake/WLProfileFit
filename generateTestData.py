import numpy as np
import NFW
import calcShapeNoise

nfw = NFW.NFW(200, 0.3, 0.7)
zl = 0.5

r_edges = np.logspace(-1., 2, 21)
r = 2./3.*(r_edges[1:]**3 - r_edges[:-1]**3)/(r_edges[1:]**2 - r_edges[:-1]**2)

M = 5
c = 4
dSigma = nfw.deltaSigma(r, M*10**14, c)

# assume errors to produce equal weight
# dSigma_err = dSigma/5.

# realistic error
dSigma_err = calcShapeNoise.calculateDSigmaErr(nfw.uni, r_edges, zl, zm = 1., n_tot = 20., sigma_e = 0.25)

np.savetxt("test.dat", np.array([r, dSigma, dSigma_err]).transpose(), header = "r dSigma dSigma_err")
