import cosmology
import numpy as np
import number_density
import NFW

def calcSigmaCr(uni, zl, chi_zs):
    sigma_cr = uni.c**2/(4*np.pi*uni.G) * cosmology.Mpc / uni.DA(zl) / (uni.Sk(chi_zs-uni.chi(zl))/uni.Sk(chi_zs)) * 10**3 / NFW.M_sun / (1.+zl)**2/10**12 # critical surface mass density [M_sun/ comoving pc^2]
    return sigma_cr

def calculateDSigmaErr(uni, l_r_edges, zl, zm = 1., n_tot = 20., sigma_e = 0.25):
    l_zs, dz = number_density.get_zs_dz()
    l_chizs = uni.chi(l_zs)

    # deltaSigma error does not depend on H0
    z0 = zm/3.
    n_s = number_density.N_gt_z(zl, z0, n_tot)

    dSigma_err = list()
    l_zs, dz = number_density.get_zs_dz()
    l_chizs = l_chizs[l_zs > zl]
    l_zs = l_zs[l_zs > zl]
    print "lens redshift:", zl
    for i in range(len(l_r_edges) - 1):
        a = 1./(1.+zl)
        area = (l_r_edges[i+1]**2 - l_r_edges[i]**2)*np.pi*(a/uni.DA(zl)[0]/np.pi*180.*60.)**2 # r is in comoving
        N = area * n_s
        l_sigma_cr = calcSigmaCr(uni, zl, l_chizs)
        l_n_s = n_tot*number_density.calcDpdz(l_zs, z0)
        l_dSigma_err = np.sqrt(sigma_e**2/(area*l_n_s*dz)*l_sigma_cr**2)
        dSigma_err_invercevariance = np.sqrt(1./(np.sum(1./l_dSigma_err**2)))
        dSigma_err.append(dSigma_err_invercevariance)

    dSigma_err = np.array(dSigma_err)
    return dSigma_err
