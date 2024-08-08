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


    #projected nfw profile
    def sigma(self, r):
        r_s = self.r_200/self.c
        k = 2*r_s*self.rho_0
        if np.isscalar(r):
           r = np.array([r])
        sig = 0.0*r
        for cnt, i in enumerate(r):
            x = i/r_s
            if x < 1:
                value = (1 - np.arccosh(1/x)/np.sqrt(1-x**2))/(x**2-1)
            elif x > 1:
                value = (1 - np.arccos(1/x)/np.sqrt(x**2-1))/(x**2-1)
            else:
                value = 1./3.
            sig[cnt] = value
        sig = value*k
        return sig



    # integrating over theta getting mean over a circle
    def avg_cir(self, r, roff):
        value = quad(lambda j: self.sigma(np.sqrt(roff**2 + r**2 + 2*roff*r*np.cos(j))), 0., 2*np.pi)[0]/(2*np.pi)
        return value

    def miscen_sigma(self, r, sigma_off):
        value= 0.0*r
        for cnt ,rr in enumerate(r):
            value[cnt] = quad(lambda j: j/sigma_off**2 * np.exp(-(j/sigma_off)**2 /2) * self.avg_cir(rr, j), 0., 10*sigma_off)[0]
        return value


if __name__ == "__main__":
    m_tot = 1e14 # m_sun
    c = 5.
    hp = halo(m_tot,4)
    rbin = np.logspace(-2,np.log10(2),10)
    plt.subplot(2,2,1)
    sigmaoff = 0.1
    plt.plot(rbin, hp.miscen_sigma(rbin, sigmaoff)/1e12,label = r'$\sigma_{\rm off} = %2.2f$'%sigmaoff)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.xlabel(r'$R [{ \rm h^{-1}Mpc}]$')
    plt.ylabel(r'$\Sigma(R) [{\rm h M_\odot pc^{-2}}]$')
    plt.savefig('offcen.png', dpi=300)


#plt.show()
#plt.plot(rd_bin,del_sigma,'r',label = 'sattelite numerical')


#c = 0
#for r in rd_bin:
#    del_sigma[c] = (2*np.pi*quad(lambda r: r*des_cir(r,roff), 0,r)[0])/(np.pi*r**2) - des_cir(r,roff)
#    print(('hello %2.5f' % (r)))
#    c = c+1
#plt.xlabel('R (Mpc)')
#plt.ylabel(r'$\Delta\Sigma (R)$')
#plt.plot(rd_bin,del_sigma,'r',label = 'sattelite numerical')
"""
kappa = np.zeros(len(rbin))
sig = np.zeros(len(rbin))

# Here I am integrated over z, r^2 = R^2 + z^2
c=0
for i  in rbin:
    ro = -np.sqrt((r_200)**2 - i**2)
    a = quad((lambda z : rho(i,z,rho_0,r_s)),ro,-ro)[0]
    kappa[c] = a
    sig[c] = sigma(i,rho_0,r_s)
    c=c+1

# interpolating the above sigma vs R.
des_pro = interp1d(rbin,kappa,kind='cubic')
del_des = np.zeros(len(rbin))
c=0
mass_in = (2*np.pi*quad(lambda x: x, 0,np.min(rbin))[0])*des_pro(np.min(rbin))

gc = 0
for i in rbin:
    del_des[gc] = (mass_in + 2*np.pi*quad(lambda x: x*des_pro(x), np.min(rbin), i)[0])/(np.pi*i**2) - des_pro(i)
    gc=gc+1


# DAUGHTER HALO PARAMETERS
# position of center of the disc
x0 = r_200/2
y0 = r_200/2
roff = np.sqrt(x0**2 + y0**2)
rd_bin = rbin/4

del_sigma = np.zeros(len(rd_bin))



    # nfw profile
    def rho(self, roff,z,rho_0,r_s):
        r = np.sqrt(roff**2 + z**2)
        des = rho_0/((r/r_s)*(1+r/r_s)**2)
        return des



"""

