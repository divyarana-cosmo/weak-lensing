import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

class constants:
    """Useful constants"""
    G = 4.301e-9 #km^2 Mpc M_sun^-1 s^-2 gravitational constant
    H0 = 100. #h kms-1 Mpc-1 hubble constant at present
    omg_m = 0.315 #omega_matter

class halo(constants):
    """Useful functions for weak lensing signal modelling"""
    def __init__(self,m_tot,con_par):
        self.m_tot = m_tot # total mass of the halo
        self.c = con_par # concentration parameter
        self.rho_crt = 3*self.H0**2/(8*np.pi*self.G) # rho critical
        self.r_200 = (3*m_tot/(4*np.pi*200*self.rho_crt*self.omg_m))**(1./3.) # radius defines size of the halo
        self.rho_0 = con_par**3 *m_tot/(4*np.pi*self.r_200**3 *(np.log(1+con_par)-con_par/(1+con_par)))
        self.init_sigma = False
        print(("Intialing NFW parameters\n m_tot = %s M_sun\nconc_parm = %s\nrho_0 = %s M_sun/Mpc^3\n r_s = %s Mpc"%(m_tot,con_par,self.rho_0,self.r_200/self.c)))

    def nfw(self,r):
        """given r, this gives nfw profile as per the instantiated parameters"""
        r_s = self.r_200/self.c
        value  = self.rho_0/((r/r_s)*(1+r/r_s)**2)
        return value

    def sigma_nfw(self,r):
        r_s = self.r_200/self.c
        k = 2*r_s*self.rho_0
        sig = 0.0*r
        c=0
        for i in r:
            x = i/r_s
            if x < 1:
                value = (1 - np.arccosh(1/x)/np.sqrt(1-x**2))/(x**2-1)
            elif x > 1:
                value = (1 - np.arccos(1/x)/np.sqrt(x**2-1))/(x**2-1)
            else:
                value = 1./3.
            sig[c] = value*k
            c=c+1

        return sig

    def avg_sigma_nfw(self,r):
        r_s = self.r_200/self.c
        k = 2*r_s*self.rho_0
        sig = 0.0*r
        c=0
        for i in r:
            x = i/r_s
            if x < 1:
                value = np.arccosh(1/x)/np.sqrt(1-x**2) + np.log(x/2.0)
                value = value*2.0/x**2
            elif x > 1:
                value = np.arccos(1/x)/np.sqrt(x**2-1)  + np.log(x/2.0)
                value = value*2.0/x**2
            else:
                value = 2*(1-np.log(2))
            sig[c] = value*k
            c=c+1

        return sig

    def esd(self,r):
        val = self.avg_sigma_nfw(r) - self.sigma_nfw(r)
        return val

    def sigma(self, R):
        if np.isscalar(R):
            Rarr = np.array([R])
        else:
            Rarr = R*1.0
        sigmaarr = Rarr*0.0
        if not self.init_sigma:
            self.init_sigma_spl()

        idx = (Rarr>=self.sigma_dict["Rmin"]) & (Rarr<self.sigma_dict["Rmax"])
        sigmaarr[idx] = 10.**self.sigma_spl(np.log10(Rarr[idx]))

        idx = (Rarr<self.sigma_dict["Rmin"])
        sigmaarr[idx] = self.sigma_dict["Sigmamin"]

        idx = (Rarr>self.sigma_dict["Rmax"])
        sigmaarr[idx] = 0.0

        if np.isscalar(R):
            sigmaarr = sigmaarr[0]

        return sigmaarr

    def init_sigma_spl(self):

        Rarr = np.logspace(-2,np.log10(8*self.r_200),30)
        Sigmaarr = Rarr*0.0

        for ii, R in enumerate(Rarr):
            #z0 = np.sqrt((8.*self.r_200)**2 - R**2)
            #if z0==0:
            #    Sigmaarr[ii]=1e-12
            #else:
            #    Sigmaarr[ii] = quad((lambda z : self.nfw(np.sqrt(R**2 + z**2))),-z0,z0)[0]
            Sigmaarr[ii] = 2*quad((lambda z : self.nfw(np.sqrt(R**2 + z**2))), 0, np.inf)[0]
        self.sigma_spl = interp1d(np.log10(Rarr), np.log10(Sigmaarr))
        self.sigma_dict = {}
        self.sigma_dict["Rmin"] = Rarr[0]
        self.sigma_dict["Rmax"] = Rarr[-1]
        self.sigma_dict["Sigmamin"] = Sigmaarr[0]
        print((Rarr[0]))
        self.init_sigma = True
        return

    def num_sigma(self,Rarr):
        if np.isscalar(Rarr):
            return 2*quad((lambda z : self.nfw(np.sqrt(Rarr**2 + z**2))), 0, np.inf)[0]

        Sigmaarr = Rarr*0.0
        for ii, R in enumerate(Rarr):
            Sigmaarr[ii] = 2*quad((lambda z : self.nfw(np.sqrt(R**2 + z**2))), 0, np.inf)[0]
        return Sigmaarr
    def num_delta_sigma(self,R):
        """difference between mean sigma and average over the circle of radius R"""
        if np.isscalar(R):
            value =  2*np.pi*quad(lambda Rp: Rp*self.num_sigma(Rp), 0.0, R)[0]/(np.pi*R**2) - self.num_sigma(R)
        else:
            value = 0.0*R
            for ii,rr in enumerate(R):
                value[ii] = 2*np.pi*quad(lambda Rp: Rp*self.num_sigma(Rp), 0.0, rr)[0]/(np.pi*rr**2) - self.num_sigma(rr)
        return value



    def delta_sigma(self,R):
        """difference between mean sigma and average over the circle of radius R"""
        if np.isscalar(R):
            value =  2*np.pi*quad(lambda Rp: Rp*self.sigma(Rp), 0.0, R)[0]/(np.pi*R**2) - self.sigma(R)
        else:
            value = 0.0*R
            for ii,rr in enumerate(R):
                value[ii] = 2*np.pi*quad(lambda Rp: Rp*self.sigma(Rp), 0.0, rr)[0]/(np.pi*rr**2) - self.sigma(rr)
        return value

    def sigma_cic(self,r,r0):
        """integrating over theta getting mean over a circle"""
        value = quad(lambda j: self.sigma(np.sqrt(r0**2 + r**2 + 2*r0*r*np.cos(j))), 0., 2*np.pi)[0]/(2*np.pi)
        return value

if __name__ == "__main__":
    plt.subplot(2,2,1)
    rbin = np.logspace(-2,np.log10(2),30)

    hp = halo(1e13,10)
    plt.plot(rbin, hp.nfw(rbin), '-', label='c=10')

    hp = halo(1e13,20)
    plt.plot(rbin, hp.nfw(rbin), '-', label='c=20')
 

    #plt.plot(rbin, hp.num_delta_sigma(rbin)/(1e12), '.', lw=0.0)
    #plt.plot(rbin, hp.esd(rbin)/(1e12))
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.xlabel(r'$r$')
    plt.ylabel(r'$\rho (r)$')
    #plt.ylabel(r'$\Delta \Sigma (R) [{\rm h M_\odot pc^{-2}}]$')
    plt.savefig('test.png', dpi=300)
