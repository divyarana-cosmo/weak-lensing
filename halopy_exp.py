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
        print("Intialing NFW parameters\n m_tot = %s M_sun\nconc_parm = %s\nrho_0 = %s M_sun/Mpc^3\n r_s = %s Mpc"%(m_tot,con_par,self.rho_0,self.r_200/self.c))

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

        print("SPLINE READY FOR NFW SIGMA")
        Rarr = np.logspace(-2,np.log10(8*self.r_200),50)
        Sigmaarr = Rarr*0.0

        for ii, R in enumerate(Rarr):
            z0 = -np.sqrt((8*self.r_200)**2 - R**2)
            if z0==0:
                Sigmaarr[ii]=self.not_so_tiny
            else:
                Sigmaarr[ii] = quad((lambda z : self.nfw(np.sqrt(R**2 + z**2))),z0,-z0)[0]
        self.sigma_spl = interp1d(np.log10(Rarr), np.log10(Sigmaarr),kind ="cubic")
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


    """segment for the parent halo contribution for the daughter halo at distance r0"""
    def sigma_cir(self,r,r0):
        """sigma mean over a circle using the spline given below"""
        #if not self.init_sigma_cir:
        #    self.init_sigma_cir_spl(r0)

        self.init_sigma_cir_spl(r0)
        if r > self.sigma_cir_dict["Rmax"]:
            value = 0.0
        elif r < self.sigma_cir_dict["Rmin"]:
            value = self.sigma_cir_dict["Sigmamin"]

        else:
            value = self.sig_cir_spl(r)

        return value

    def init_sigma_cir_spl(self,r0):
        """spline for the satellite at a distance r0 from the center for parents contribution averaged over a circle"""

        print("SPLINE READY FOR AVERAGING OVER CIRCLE")
        rdbin = np.logspace(-2,np.log(10*self.r_200),50)
        des_cir = 0.0*rdbin
        for i  in range(0,len(rdbin)):
            des_cir[i] = quad(lambda j: self.sigma(np.sqrt(r0**2 + rdbin[i]**2 + 2*r0*rdbin[i]*np.cos(j))), 0., 2*np.pi)[0]/(2*np.pi)

        self.sig_cir_spl = interp1d(rdbin,des_cir,kind = "cubic")
        self.sigma_cir_dict = {}
        self.sigma_cir_dict["Rmax"] = rdbin[-1]
        self.sigma_cir_dict["Rmin"] = rdbin[0]
        self.sigma_cir_dict["Sigmamin"] = des_cir[0]
        self.init_sigma_cir = True
        return

    def delta_sigma_dau(self,r,r0):
        value =  2*np.pi*quad(lambda rp: rp*self.sigma_cir(rp,r0), 0.0, r)[0]/(np.pi*r**2) - self.sigma_cir(r,r0)
        return value



if __name__ == "__main__":
    """mtot_bin = np.arange(10,11,1)
    mtot_bin = 10**mtot_bin
    #c_bin = np.
    rdbin = np.logspace(-1,np.log10(10),50)
    for m in mtot_bin:
        h_p = halo(m,10.0)
        delta_part = 0.0*rdbin
        for i in range(0,len(rdbin)):
            if rdbin[i]>=(8*h_p.r_200):
                delta_part[i] = h_p.delta_sigma(8*h_p.r_200) + h_p.sigma(8*h_p.r_200)
            else:
                delta_part[i] = h_p.delta_sigma(rdbin[i])

        #print rdbin[i]

"""
    #print h_d.delta_sigma(rdbin[0])
    #print h_p.delta_sigma_dau(rdbin[0],rd_dist)
    hp = halo(1e12,10.0)
    r0 = hp.r_200/2
    rbin = np.logspace(-2,np.log10(4),10)
    dcic = rbin*0.0
    for i in range(0,len(rbin)):
        dcic[i] = hp.sigma_cir(rbin[i],r0)

    plt.plot(rbin,dcic,'or')

    r0 = hp.r_200/4
    dcic1 = rbin*0.0
    for i in range(0,len(rbin)):
        dcic1[i] = hp.sigma_cir(rbin[i],r0)

    plt.plot(rbin,dcic1,label='sdhja')
    #print delta_part
        #delta_part = 0.0*rdbin

    #plt.plot(rdbin,delta_part/1e12,'ob', label = 'parent')
    #plt.plot(rdbin,(delta_daug + delta_part)/1e12,'og', label  = 'addition of both')
    #plt.axvline(6*h_p.r_200)
    #plt.axvline(2.0*rd_dist)
    #plt.plot(rdbin,h_p.sigma(rdbin))
    #plt.xlim(0.1,4)
    #plt.ylim(1e0,)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$R (Mpc h^{-1})$')
    plt.ylabel(r'$\Delta \Sigma (R) \times 10^{12} M_\odot (Mpc/h)^{-2}$ ')
    #plt.savefig("/Users/divyarana/Dropbox/grad_school_project/images/del_sigma.pdf",format='pdf')
    plt.legend()
    plt.show()
