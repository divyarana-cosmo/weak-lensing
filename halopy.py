import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
class constants:
    """Useful constants"""
    G = 4.301e-9 #km^2 Mpc M_sun^-1 s^-2 gravitational constant
    H0 = 67. #kms-1 Mpc-1 hubble constant at present
    omg_m = 0.3 #omega_matter

class halo(constants):
    """Useful functions for weak lensing signal modelling"""
    def __init__(self,m_tot,con_par):
        self.m_tot = m_tot # total mass of the halo
        self.c = con_par # concentration parameter
        self.rho_crt=3*self.H0**2/(8*np.pi*self.G)
        self.r_200 = (3*m_tot/(4*np.pi*200*self.rho_crt*self.omg_m))**(1./3.)
        self.rho_0 = con_par**3 *m_tot/(4*np.pi*self.r_200**3 *(np.log(1+con_par)-con_par/(1+con_par)))

        print "Intialing NFW parameters\n m_tot = %s M_sun\nconc_parm = %s\nrho_0 = %s M_sun/Mpc^3\n r_s = %s Mpc"%(m_tot,con_par,self.rho_0,self.r_200/self.c)

    def nfw(self,r):
        """given r, this gives nfw profile as per the instantiated parameters"""
        r_s = self.r_200/self.c
        value  = self.rho_0/((r/r_s)*(1+r/r_s)**2)
        return value

    def sigma(self,R):
        """Projected density at a distance R from the center"""
        if(R <=0.01):
            R = 0.01
            z0 = -np.sqrt(self.r_200**2 - R**2)
            value = quad((lambda z : self.nfw(np.sqrt(R**2 + z**2))),z0,-z0)[0]
        else:
            z0 = -np.sqrt(self.r_200**2 - R**2)
            value = quad((lambda z : self.nfw(np.sqrt(R**2 + z**2))),z0,-z0)[0]
        return value

    def delta_sigma(self,R):
        """difference between mean sigma and average over the circle of radius R"""
        rbin = np.linspace(0.0,self.r_200,100)
        sig = np.zeros(len(rbin))
        for i  in range(0,len(rbin)):
            sig[i] = self.sigma(rbin[i])
        sig_of_r = interp1d(rbin,sig,kind='cubic')

        value =  2*np.pi*quad(lambda r: r*self.sigma(r), 0.0, R)[0]/(np.pi*R**2) - self.sigma(R)
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
