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
        self.rho_0 = con_par**3 *m_tot/(4*np.pi*self.r_200**3 *(np.log(1.0+con_par)-con_par/(1.0+con_par)))
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

    def sigma_cic(self,r,r0):
        """integrating over theta getting mean over a circle"""
        value = quad(lambda j: self.sigma(np.sqrt(r0**2 + r**2 + 2*r0*r*np.cos(j))), 0., 2*np.pi)[0]/(2*np.pi)
        return value

if __name__ == "__main__":
    """initiating the class"""
    #contribution due to daughter nfw
    h_d = halo(2e12,10.)


    rdbin = np.logspace(-2.6,np.log10(4*h_d.r_200),50)
    delta_sig_dau_nfw = 0.0*rdbin
    for i  in range(0,len(rdbin)):
        delta_sig_dau_nfw[i] = h_d.delta_sigma(rdbin[i])

    #contribution due to parent halo nfw
    h_p = halo(2e14,10.)
    print "here"
    print h_p.sigma(0.4)

    """position of the center of the daughter halo"""
    x0 = h_p.r_200/2
    y0 = h_p.r_200/2
    r0 = np.sqrt(x0**2 + y0**2)
    #rbin = np.logspace(-2,np.log10(h_p.r_200),20)
    des_cir = 0.0*rdbin
    for i  in range(0,len(rdbin)):
        des_cir[i] = h_p.sigma_cic(rdbin[i],r0)

    sig_cir_spl = interp1d(rdbin,des_cir,kind = "cubic")
    def sig_cir(r):
        if r < np.min(rdbin):
            value =  sig_cir_spl(np.min(rdbin))
        elif r > np.max(rdbin):
            value = 0
        else:
            value = sig_cir_spl(r)
        return value
    dau_delta_sigma = 0.0*rdbin
    c = 0
    for r1 in rdbin:
        dau_delta_sigma[c] = (2*np.pi*quad(lambda r: r*sig_cir(r), 0,r1)[0])/(np.pi*r1**2) - sig_cir(r1)
        print 'hello %2.5f' % (r1)
        c = c+1


    #plt.plot(rdbin,dau_delta_sigma/1e12 ,'r', label = 'parent halo contribution')
    #plt.plot(rdbin,delta_sig_dau_nfw/1e12,'b', label  = 'daughter halo contribution')
    #plt.plot(rdbin,(delta_sig_dau_nfw + dau_delta_sigma)/1e12,'g', label  = 'addition of both')
    #plt.axvline(r0)
    plt.plot(rdbin,h_p.sigma(rdbin))
    #plt.xlim(0.05,)
    plt.xscale('log')
    plt.xlabel('R (Mpc h-1)')
    plt.ylabel(r'$\Delta \Sigma (R) \times 10^{12}$  ')

    plt.yscale('log')
    plt.legend()
    plt.show()
