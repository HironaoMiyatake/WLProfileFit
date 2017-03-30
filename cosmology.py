import numpy
from scipy import integrate
from scipy import interpolate

Mpc = 3.08568025*10**19 # 1Mpc = 3.08568025*10**19 km
M_sun = 1.9884*10**33 # solar mass [g] 

class universe:
    splineOneOverH = None
    
    def __init__(self, Omega_m0 = 0.3, Omega_l0 = 0.7, K = 0, H0 = 70., c = 2.99792458*10**5, G = 6.67384*10**(-20), z_step = 0.000001, z_max = 10.):
        self.Omega_m0 = Omega_m0
        self.Omega_l0 = Omega_l0
        self.K = K
        self.H0 = H0 # Hubble paramter at present [km/s/Mpc]
        self.c = c # speed of light [km/s]
        self.G = G # gravitational constant [km^3s^-2kg^-1]
        self.rho_cr0 = (3.*H0**2)/(8.*numpy.pi*G)/Mpc**2 * 10**-12 # critical density at present [g/cm^3]
        self.z_step = z_step # step for integration relavant to z
        self.z_max = z_max

    # Hubble parameter
    def H(self, z):
        return self.H0*numpy.sqrt(self.Omega_m0*numpy.power(1+z,3) + self.K*self.c**2/self.H0/self.H0*numpy.power(1+z,2) + self.Omega_l0)

    def H_a(self, a):
        return self.H0*numpy.sqrt(self.Omega_m0/numpy.power(a, 3) + self.K*self.c**2/self.H0/self.H0/numpy.power(a,2) + self.Omega_l0)

    # E(z) = H(z)/H0
    def E(self, z):
        return self.H(z)/self.H0

    def E_a(self, a):
        return self.H_a(a)/self.H0

    def OneOverH(self, z):
        return 1./self.H(z)

    def calcSplineOneOverH(self):
        zz_crude = numpy.linspace(0., self.z_max, 10000.)
        oneOverH_crude = self.OneOverH(zz_crude)
        self.splineOneOverH = interpolate.splrep(zz_crude, oneOverH_crude, s=0)

    # comoving distance
    def chi_old(self, z, n_bin = 1000):
        def integrateOneOverH(z, n_bin = 1000):
            zz = numpy.linspace(0., z, n_bin)
            interval = zz[1] - zz[0]
            return self.OneOverH(zz[0])*0.5*interval + self.OneOverH(zz[1:n_bin-1]).sum()*interval + self.OneOverH(zz[n_bin-1])*0.5*interval
            
        if isinstance(z,list) or isinstance(z,numpy.ndarray):
            chi = list()
            i = 0
            for zz in z:
                chi.append(self.c*integrateOneOverH(zz, n_bin))
            return numpy.array(chi)
        else:
            chi = self.c*integrateOneOverH(z)
            return chi

    def chi(self, z, n_bin = 1000):
        if (isinstance(z, list) or isinstance(z,numpy.ndarray)) == False:
            z = numpy.array([z])
        zzz = numpy.empty((len(z), n_bin))
        interval = list()
        for i, zz in enumerate(z):
            z_tmp, step = numpy.linspace(0., zz, n_bin, retstep = True)
            zzz[i] = z_tmp
            interval.append(step)
        interval = numpy.array(interval)
        oneOverH = self.OneOverH(zzz)
        return self.c*((oneOverH[:,0] + oneOverH[:,n_bin-1])*0.5 + oneOverH[:,1:n_bin-1].sum(axis = 1))*interval

    try:
        from scipy import integrate
        from scipy import interpolate
    except ImportError:
        pass
    else:
        def chi_scipy(self, z):
            if isinstance(z,list) or isinstance(z,numpy.ndarray):
                chi = list()
                i = 0
                for zz in z:
                    y, abserr = integrate.quad(self.OneOverH, 0., zz)
                    chi.append(y*self.c)
                return numpy.array(chi)
            else:
                chi = self.c*integrate.quad(self.OneOverH, 0., z)[0]
                return chi

        def chi_scipy_2(self, z, n_bin = 1000):
            def integrateOneOverH(z, n_bin = 1000):
                if z > self.z_max:
                    zz_crude = numpy.linspace(0., z, 10000.)
                    oneOverH_crude = self.OneOverH(zz_crude)
                    sp = interpolate.splrep(zz_crude, oneOverH_crude, s=0)
                else:
                    if not self.splineOneOverH:
                        self.calcSplineOneOverH()
                    sp = self.splineOneOverH

                zz = numpy.linspace(0., z, n_bin)
                oneOverH = interpolate.splev(zz, sp, der=0)
                interval = zz[1] - zz[0]
                return oneOverH[0]*0.5*interval + oneOverH[1:n_bin-1].sum()*interval + oneOverH[n_bin-1]*0.5*interval

            if isinstance(z,list) or isinstance(z,numpy.ndarray):
                chi = list()
                i = 0
                for zz in z:
                    chi.append(self.c*integrateOneOverH(zz, n_bin))
                return numpy.array(chi)
            else:
                chi = self.c*integrateOneOverH(z, n_bin)
                return chi

        # Growing solution of matter over density delta(a) \prop D
        def D_a(self, a): # it is not normalized
            def integrand_D_a(a):
                return 1./a**3/self.H_a(a)**3*self.H0**3
            return 5./2.*self.Omega_m0*self.H_a(a)/self.H0*integrate.quad(integrand_D_a, 0., a)[0]

        # Growth rate f(z) = dlnD/dlna
        def f_a(self, a): # it is normilized by definition
            return -0.5*(self.Omega_m0*a**-3-2.*self.Omega_l0)/self.E_a(a)**2 - 1. + 2.5*self.Omega_m0*a**-2/self.E_a(a)**2/self.D_a(a)
            #return -3.*a**-3*self.Omega_m0/2./self.E_a(a)**2 + 2.5*self.Omega_m0*a**-2/self.E_a(a)**2/self.D_a(a)
        

    # Sk
    def Sk(self, chi):
        if self.K > 0:
            return 1./numpy.sqrt(self.K)* numpy.sin(chi*numpy.sqrt(self.K))
        elif self.K == 0:
            return chi
        elif self.K < 0:
            return 1./numpy.sqrt(-self.K)* numpy.sinh(chi*numpy.sqrt(-self.K))

    # angular diameter distance
    def DA(self, z):
        return self.Sk(self.chi(z))/(1+z)
        
    # mean mass density at redshift z
    def rho_m(self, z):
        return self.Omega_m0 * self.rho_cr0 * numpy.power(1+z,3)

    # critical density at redshift z
    def rho_cr(self, z):
        return (3.*self.H(z)**2)/(8.*numpy.pi*self.G)/Mpc**2 * 10**-12
