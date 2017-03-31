import numpy as np
import matplotlib.pyplot as plt

# dndz from Eq. 17 in Oguri & Takada (2011)

def calcDpdz(zs, z0):
    return zs**2/2./z0**3*np.exp(-zs/z0)

def N_gt_z(zs_min, z0, n_tot):
    zs, dz = np.linspace(0., 5., 10000, retstep = True)
    return n_tot*dz*np.sum(calcDpdz(zs, z0)[zs > zs_min])

def zs_gt_z(zs_min, z0):
    zs, dz = np.linspace(0., 5., 10000, retstep = True)
    return np.sum(zs[zs > zs_min]*calcDpdz(zs, z0)[zs > zs_min])/np.sum(calcDpdz(zs, z0)[zs > zs_min])

def get_zs_dz():
    zs, dz = np.linspace(0., 5., 10000, retstep = True)
    return zs, dz

