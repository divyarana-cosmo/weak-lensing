import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
class constants:
    """Useful constants"""
    G = 4.301e-9 #km^2 Mpc M_sun^-1 s^-2 gravitational constant
    H0 = 100. #kms-1 Mpc-1 hubble constant at present
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
        print "Intialing NFW parameters\n m_tot = %s M_sun\nconc_parm = %s\nrho_0 = %s M_sun/Mpc^3\n r_s = %s Mpc"%(m_tot,con_par,self.rho_0,self.r_200/self.c)

    def nfw(self,r):
        """given r, this gives nfw profile as per the instantiated parameters"""
        r_s = self.r_200/self.c
        value  = self.rho_0/((r/r_s)*(1+r/r_s)**2)
        return value

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

        Rarr = np.logspace(-2,np.log10(5*self.r_200),30)
        Sigmaarr = Rarr*0.0

        for ii, R in enumerate(Rarr):
            z0 = -np.sqrt((8.*self.r_200)**2 - R**2)
            Sigmaarr[ii] = quad((lambda z : self.nfw(np.sqrt(R**2 + z**2))),z0,-z0)[0]
        self.sigma_spl = interp1d(np.log10(Rarr), np.log10(Sigmaarr))
        self.sigma_dict = {}
        self.sigma_dict["Rmin"] = Rarr[0]
        self.sigma_dict["Rmax"] = Rarr[-1]
        self.sigma_dict["Sigmamin"] = Sigmaarr[0]
        self.init_sigma = True
        return

    def delta_sigma(self,R):
        """difference between mean sigma and average over the circle of radius R"""

        value =  2*np.pi*quad(lambda Rp: Rp*self.sigma(Rp), 0.0, R)[0]/(np.pi*R**2) - self.sigma(R)

        return value


if __name__ == "__main__":
    h_p = halo(2e14,10.)

    rbin = np.logspace(-2,np.log10(h_p.r_200),20)
    sig = np.zeros(len(rbin))
    for i  in range(0,len(rbin)):
        sig[i] = h_p.delta_sigma(rbin[i])

    plt.plot(rbin,sig)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
