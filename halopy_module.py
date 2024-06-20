import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

class constants:
    """Useful constants"""
    G = 4.301e-9 #km^2 Mpc M_sun^-1 s^-2 gravitational constant
    H0 = 100. #h kms-1 Mpc-1 hubble constant at present
    omg_m = 0.315 #omega_matter
    not_so_tiny = 1e-24
class halo(constants):
    """Useful functions for weak lensing signal modelling"""
    def __init__(self,m_tot,con_par):
        self.m_tot = m_tot # total mass of the halo
        self.c = con_par # concentration parameter
        self.rho_crt = 3*self.H0**2/(8*np.pi*self.G) # rho critical
        self.r_200 = (3*m_tot/(4*np.pi*200*self.rho_crt*self.omg_m))**(1./3.) # radius defines size of the halo
        self.rho_0 = con_par**3 *m_tot/(4*np.pi*self.r_200**3 *(np.log(1.0+con_par)-con_par/(1.0+con_par)))
        self.init_sigma = False
        self.init_sigma_cir = False
        self.sigma_cir_dict = {}
        print(("Intialing NFW parameters\n m_tot = %s M_sun\nconc_parm = %s\nrho_0 = %s M_sun/Mpc^3\n r_s = %s Mpc"%(m_tot,con_par,self.rho_0,self.r_200/self.c)))

    def nfw(self,r):
        """given r, this gives nfw profile as per the instantiated parameters"""
        r_s = self.r_200/self.c
        value  = self.rho_0/((r/r_s)*(1+r/r_s)**2)
        return value

    def sigma(self, r):
        r_s = self.r_200/self.c
        k = 2*r_s*self.rho_0
        if np.isscalar(r):
           r = np.array([r])
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
        c=0
        if np.isscalar(r):
            r = np.array([r])
        sig = 0.0*r
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


    def delta_sigma(self,r):
        """difference between mean sigma and average over the circle of radius R"""
        val = self.avg_sigma_nfw(r) - self.sigma(r)
        return val


    """segment for the parent halo contribution for the daughter halo at distance r0"""
    def sigma_cir(self,r,r0):
        """sigma mean over a circle using the spline given below"""
        if not self.init_sigma_cir:
            self.init_sigma_cir_spl(r0)

        if r > self.sigma_cir_dict["Rmax"]:
            value = quad(lambda j: self.sigma(np.sqrt(r0**2 + r**2 + 2*r0*r*np.cos(j))), 0., 2*np.pi)[0]/(2*np.pi)
        elif r < self.sigma_cir_dict["Rmin"]:
            value = self.sigma_cir_dict["Sigmamin"]
        else:
            value = 10**self.sig_cir_spl(np.log10(r))
        return value

    def init_sigma_cir_spl(self,r0):
        """spline for the satellite at a distance r0 from the center for parents contribution averaged over a circle"""

        print("SPLINE READY FOR AVERAGING OVER CIRCLE")
        rdbin = np.logspace(-3,np.log10(10*self.r_200),50)
        des_cir = 0.0*rdbin
        for i  in range(0,len(rdbin)):
            des_cir[i] = quad(lambda j: self.sigma(np.sqrt(r0**2 + rdbin[i]**2 + 2*r0*rdbin[i]*np.cos(j))), 0., 2*np.pi)[0]/(2*np.pi)

        self.sig_cir_spl = interp1d(np.log10(rdbin), np.log10(des_cir),kind = "cubic")
        self.sigma_cir_dict["Rmax"] = rdbin[-1]
        self.sigma_cir_dict["Rmin"] = rdbin[0]
        self.sigma_cir_dict["Sigmamin"] = des_cir[0]
        self.init_sigma_cir = True
        return

    def delta_sigma_dau(self,r,r0):
        value =  2*np.pi*quad(lambda rp: rp*self.sigma_cir(rp,r0), 0.0, r)[0]/(np.pi*r**2) - self.sigma_cir(r,r0)
        #value =  self.sigma_cir(r,r0)
        return value



if __name__ == "__main__":
    rdbin = np.logspace(-2,np.log10(5),50)
    mhpart = 1e14
    mhdaut = 1e12
    msteldaut = 1e10

    h_p = halo(mhpart,4)
    h_d = halo(mhdaut,4)
    rd_dist =  0.3
    delta_part = 0.0*rdbin
    for i in range(len(rdbin)):
        delta_part[i] = h_p.delta_sigma_dau(rdbin[i], rd_dist)

    plt.subplot(2,2,1)
    plt.plot(rdbin,delta_part/1e12,'-', label = 'parent')
    plt.plot(rdbin,(h_d.delta_sigma(rdbin) + msteldaut/(np.pi*rdbin**2))/1e12,'-', label = 'daughter')
    plt.plot(rdbin,(delta_part + msteldaut/(np.pi*rdbin**2) + h_d.delta_sigma(rdbin))/1e12,'-', label = 'total')
    plt.xscale('log')
    plt.legend()
    plt.axhline(0.0, ls='--',color='grey')
    plt.axvline(rd_dist, color='black')
    plt.xlabel(r'$R [{ \rm h^{-1}Mpc}]$')
    plt.ylabel(r'$\Delta \Sigma [{\rm h M_\odot pc^{-2}}]$')
    plt.savefig('test.png', dpi=300)

