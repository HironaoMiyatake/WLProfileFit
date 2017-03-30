import cosmology
import scipy.optimize as spOptimize
import scipy.special as spSpecial
import numpy as np

M_sun = 1.9884*10**33 # solar mass [g] 

# Mass is M_delta_rho_m0 = M_sun/h.
# input r in units of comoving Mpc/h

def convertC(c, delta, target_delta):
    def F(target_c):
        return target_delta*target_c**3 - delta*c**3*(np.log(1.+c) - c/(1. + c))**-1 * (np.log(1.+target_c) - target_c/(1.+target_c))
    target_c = spOptimize.newton(F, c)
    return target_c

class NFW:
    def __init__(self, delta, Omega_m0, Omega_l0, cMRelation = None):
        self.delta = delta
        self.uni = cosmology.universe(Omega_m0 = Omega_m0,
                                           Omega_l0 = Omega_l0,
                                           H0 = 100.
                                           )
        self.mass_base = self.uni.rho_m(0.)/M_sun * (cosmology.Mpc*10**5)**3 
        if cMRelation is None:
            self.cMRelation = None
        else:
            if cMRelation[0] == "RM04":
                self.cMRelation = cMRelation[0]
                import sys
                sys.path.append("../model")
                import extendedCosmology
                import Mandelbaum04 as M04
                root = cMRelation[1]
                sigma8 = cMRelation[2]
                z = cMRelation[3]
                e = extendedCosmology.extendedCosmology("../model/%s_params.ini" % root, "../model/%s_matterpower.dat" % root, sigma8 = sigma8)
                self.Mnl = M04.Mnl(e, z)
                print "Use c-M relation in RM04, Mnl: %s[Msun/h]" % self.Mnl
            elif cMRelation[0] == "RM08":
                self.cMRelation = cMRelation[0]
                self.z = cMRelation[1]
                print "Use c-M relation in RM08"

    def r_s(self, M_delta, c_delta):
        return np.power(3*M_delta/(4*np.pi*self.mass_base*self.delta),1./3.)/c_delta

    def rho_s(self, M_delta, c_delta):
        return c_delta**3*self.mass_base*self.delta/3./(np.log(1.+c_delta)-c_delta/(1.+c_delta))

    def rho(self, r, M_delta, c_delta):
        rho_s = self.rho_s(M_delta, c_delta)
        r_s = self.r_s(M_delta, c_delta)
        rho = rho_s/(r/r_s)/(1.+r/r_s)**2
        return rho

    def BMO_rho(self, r, M_delta, c_delta, tau_v=2.6):
        rho_s = self.rho_s(M_delta, c_delta)
        r_s = self.r_s(M_delta, c_delta)
        rt = self.r_s(M_delta, c_delta)*c_delta*tau_v
        rho = rho_s/(r/r_s)/(1.+r/r_s)**2*(rt**2/(r**2+rt**2))**2
        return rho

    # lensing profile
    def deltaSigma(self, r, M_delta, c_delta):
        if self.cMRelation == None:
            _c_delta = c_delta
        elif self.cMRelation == "RM04":
            _c_delta = self.calcCFromCMrelationRM04(M_delta)
        elif self.cMRelation == "RM08":
            _c_delta = self.calcCFromCMrelationRM08(M_delta)

        r_s = self.r_s(M_delta, _c_delta)

        x = r/r_s
        x_ltone = x[x < 1]
        x_eqone = x[np.equal(x, 1.)]
        x_gtone = x[x > 1]
        
        fact = r_s*self.rho_s(M_delta, _c_delta)

        deltaSigma_ltone = fact * (8.*np.arctanh(np.sqrt((1.-x_ltone)/(1.+x_ltone)))/(x_ltone**2*np.sqrt(1.-x_ltone**2))
                              +4.*np.log(x_ltone/2.)/x_ltone**2
                              -2./(x_ltone**2-1.)
                              +4.*np.arctanh(np.sqrt((1.-x_ltone)/(1.+x_ltone)))/((x_ltone**2-1.)*np.sqrt(1.-x_ltone**2)))
        deltaSigma_eqone = fact * (10./3.+4.*np.log(0.5)) * x_eqone
        deltaSigma_gtone = fact * (8.*np.arctan(np.sqrt((x_gtone-1.)/(1.+x_gtone)))/x_gtone**2/np.sqrt(x_gtone**2-1.)
                              +4.*np.log(x_gtone/2.)/x_gtone**2
                              -2./(x_gtone**2-1.)
                              +4.*np.arctan(np.sqrt((x_gtone-1.)/(1.+x_gtone)))/(x_gtone**2-1.)**1.5)
        deltaSigma = np.concatenate((deltaSigma_ltone, deltaSigma_eqone, deltaSigma_gtone))/10**12 # 10**12 converts [hM_sun/Mpc^2 to hM_sun/pc^2]

        return deltaSigma

    def truncated_deltaSigma(self, r, M_delta, c_delta):
        r_s = self.r_s(M_delta, c_delta)
        r_delta = r_s*c_delta

        x = r/r_s

        x_ltone = x[x < 1]
        x_eqone = x[np.equal(x, 1.)]
        x_gtoneltc = x[np.logical_and(x > 1, x <= c_delta)]
        x_gtc = x[x > c_delta]

        fact = M_delta/(np.log(1.+c_delta)-c_delta/(1.+c_delta))*c_delta**2/2./np.pi/r_delta**2

        deltaSigma_ltone = 1./(x_ltone**2*(1.+c_delta))*((2.-x_ltone**2)*np.sqrt(c_delta**2-x_ltone**2)/(1-x_ltone**2)-2.*c_delta) + 2/x_ltone**2*np.log(x_ltone*(1.+c_delta)/(c_delta+np.sqrt(c_delta**2-x_ltone**2))) + (2.-3*x_ltone**2)/(x_ltone**2*(1.-x_ltone**2)**1.5)*np.arccosh((x_ltone**2+c_delta)/(x_ltone*(1.+c_delta)))
        deltaSigma_eqone = 1./3./(1.+c_delta)*((11*c_delta+10)*np.sqrt(c_delta**2-1.)/(1.+c_delta)-6*c_delta) + 2.*np.log((1.+c_delta)/(c_delta+np.sqrt(c_delta**2-1.))) * x_eqone
        deltaSigma_gtoneltc = 1./(x_gtoneltc**2*(1.+c_delta))*((2.-x_gtoneltc**2)*np.sqrt(c_delta**2-x_gtoneltc**2)/(1-x_gtoneltc**2)-2.*c_delta) + 2/x_gtoneltc**2*np.log(x_gtoneltc*(1.+c_delta)/(c_delta+np.sqrt(c_delta**2-x_gtoneltc**2))) - (2.-3*x_gtoneltc**2)/(x_gtoneltc**2*(x_gtoneltc**2-1.)**1.5)*np.arccos((x_gtoneltc**2+c_delta)/(x_gtoneltc*(1.+c_delta)))
        deltaSigma_gtc = 2/x_gtc**2*(np.log(1.+c_delta)-c_delta/(1.+c_delta))
        deltaSigma = fact*np.concatenate((deltaSigma_ltone, deltaSigma_eqone, deltaSigma_gtoneltc, deltaSigma_gtc))/10**12 # 10**12 converts [hM_sun/Mpc^2 to hM_sun/pc^2]
        return deltaSigma

    def BMO_deltaSigma(self, r, M_delta, c_delta, tau_v=2.6):
        rho_s = self.rho_s(M_delta, c_delta)
        r_s = self.r_s(M_delta, c_delta)

        tau = tau_v*c_delta
        x = r/r_s

        x_ltone = x[x < 1.]
        F_ltone = 1./np.sqrt(1.-x_ltone**2)*np.arctanh(np.sqrt(1.-x_ltone**2))
        x_gtone = x[x > 1.]
        F_gtone = 1./np.sqrt(x_gtone**2-1.)*np.arctan(np.sqrt(x_gtone**2-1.))
        F = np.concatenate([F_ltone, F_gtone])

        L = np.log(x/(np.sqrt(tau**2+x**2)+tau))

        kappa_sigcrit = rho_s*r_s*tau**4/(tau**2+1.)**3*(2.*(tau**2+1.)/(x**2-1)*(1.-F)+8.*F+(tau**4-1.)/tau**2/(tau**2+x**2)-np.pi*(4.*(tau**2+x**2)+tau**2+1)/(tau**2+x**2)**1.5+(tau**2*(tau**4-1.)+(tau**2+x**2)*(3.*tau**4-6.*tau**2-1.))/tau**3/(tau**2+x**2)**1.5*L)
        kappabar_sigcrit = 2.*rho_s*r_s*tau**4/(tau**2+1.)**3/x**2*(2.*(tau**2+4.*x**2-3.)*F+1./tau*(np.pi*(3.*tau**2-1.)+2.*tau*(tau**2-3.)*np.log(tau))+1./tau**3/np.sqrt(tau**2+x**2)*(-tau**3*np.pi*(4.*x**2+3.*tau**2-1.)+(2.*tau**4*(tau**2-3.)+x**2*(3.*tau**4-6.*tau**2-1.))*L))

        deltaSigma = (kappabar_sigcrit - kappa_sigcrit)/10**12
        
        return deltaSigma

    def BMO_deltaSigma_kappacorr(self, r, M_delta, c_delta, Sigma_cr, tau_v=2.6):
        rho_s = self.rho_s(M_delta, c_delta)
        r_s = self.r_s(M_delta, c_delta)

        tau = tau_v*c_delta
        x = r/r_s

        x_ltone = x[x < 1.]
        F_ltone = 1./np.sqrt(1.-x_ltone**2)*np.arctanh(np.sqrt(1.-x_ltone**2))
        x_gtone = x[x > 1.]
        F_gtone = 1./np.sqrt(x_gtone**2-1.)*np.arctan(np.sqrt(x_gtone**2-1.))
        F = np.concatenate([F_ltone, F_gtone])

        L = np.log(x/(np.sqrt(tau**2+x**2)+tau))

        kappa_sigcrit = rho_s*r_s*tau**4/(tau**2+1.)**3*(2.*(tau**2+1.)/(x**2-1)*(1.-F)+8.*F+(tau**4-1.)/tau**2/(tau**2+x**2)-np.pi*(4.*(tau**2+x**2)+tau**2+1)/(tau**2+x**2)**1.5+(tau**2*(tau**4-1.)+(tau**2+x**2)*(3.*tau**4-6.*tau**2-1.))/tau**3/(tau**2+x**2)**1.5*L)
        kappabar_sigcrit = 2.*rho_s*r_s*tau**4/(tau**2+1.)**3/x**2*(2.*(tau**2+4.*x**2-3.)*F+1./tau*(np.pi*(3.*tau**2-1.)+2.*tau*(tau**2-3.)*np.log(tau))+1./tau**3/np.sqrt(tau**2+x**2)*(-tau**3*np.pi*(4.*x**2+3.*tau**2-1.)+(2.*tau**4*(tau**2-3.)+x**2*(3.*tau**4-6.*tau**2-1.))*L))

        deltaSigma = (kappabar_sigcrit - kappa_sigcrit)/10**12
        kappa = kappa_sigcrit/10**12/Sigma_cr
        deltaSigma_corr = deltaSigma*(1+kappa)
        
        return deltaSigma_corr

    def BMO_uk(self, k, M_delta, c_delta, tau_v=2.6):
        r_s = self.r_s(M_delta, c_delta)
        tau = tau_v*c_delta
        x = k*r_s
        mnfw = np.log(1.+c_delta)-c_delta/(1.+c_delta)

#        shi_taux, chi_taux = spSpecial.shichi(tau*x)
#        P_taux = np.sinh(tau*x)*chi_taux-np.cosh(tau*x)*shi_taux
#        Q_taux = np.cosh(tau*x)*chi_taux-np.sinh(tau*x)*shi_taux

        gamma = 0.57721566
        def Pfit(x):
            a,b,c,d,e,f = 1.5652, 3.38723, 6.34891, 0.817677, -0.0895584, 0.877375
            return -(1./x+(b*x**e/(c+(x-d)**2)))*(x**4/(x**4+a**4))**f+x*(gamma+np.log(x)-1.)*(a**4/(x**4+a**4))**f
        def Qfit(x):
            a,b,c,d,e,f,g = 2.26901, -2839.04, 265.511, -1.12459, -2.90136, 1.86475, 1.52197
            return (1./x**2+(b*x**e/(c+(x-d)**4)))*(x**4/(x**4+a**4))**g+((gamma+np.log(x))*(1.+0.5*x**2)-3./4.*x**2)*(a**4/(x**4+a**4))**f

        si_x, ci_x = spSpecial.sici(x)
        uk = tau/4./mnfw/(1+tau**2)**3/x*(2.*(3.*tau**4-6.*tau**2-1.)*Pfit(tau*x)-2.*tau*(tau**4-1.)*x*Qfit(tau*x)-2.*tau**2*np.pi*np.exp(-tau*x)*((tau**2+1)*x+4.*tau)+2.*tau**3*(np.pi-2.*si_x)*(4.*np.cos(x)+(tau**2+1)*x*np.sin(x))+4.*tau**3*ci_x*(4.*np.sin(x)-(tau**2+1)*x*np.cos(x)))

        return uk

    def calcCFromCMrelationRM04(self, M_delta):
        if self.delta != 180.:
            raise ValueError, "delta should be 180 to use CM relation"
        return 10.*(M_delta/self.Mnl)**-0.13

    def calcCFromCMrelationRM08(self, M_delta):
        delta200 = 200.
        M0 = 1.e14
        f = 1.22/(1.+self.z)*4.6
        beta = -0.13
        if self.delta != 200.:
            G = M0/delta200/f**(1./beta)
            def F(c):
                c200 = (M_delta/c**3/self.delta/G)**(beta/(1.-3*beta))
                return self.delta*c**3 - delta200*c200**3*(np.log(1.+c200)-c200/(1.+c200))**-1*(np.log(1.+c)-c/(1.+c)) 
            c_delta = spOptimize.newton(F, 6.)
            c200 = (M_delta/c_delta**3/self.delta/G)**(beta/(1.-3*beta))
            M200 = (c200/f)**(1./beta) * M0
            #print "M200/c200**3/delta200:", M200/c200**3/delta200
            #print "M_delta/c_delta**3/delta:", M_delta/c_delta**3/self.delta
            return c_delta
        else:
            c_delta = f*(M_delta/M0)**beta
            return c_delta

