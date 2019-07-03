import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d

# nfw profile
def rho(r0,z,rho_0,r_s):
    r = np.sqrt(r0**2 + z**2)
    des = rho_0/((r/r_s)*(1+r/r_s)**2)
    return des

#projected nfw profile
def sigma(r,rho_0,r_s):
    k = 2*r_s*rho_0
    x = r/r_s
    if x < 1:
        value = (1 - np.arccosh(1/x)/np.sqrt(1-x**2))/(x**2-1)
    elif x > 1:
        value = (1 - np.arccos(1/x)/np.sqrt(x**2-1))/(x**2-1)
    else:
        value = 1./3.

    value = value*k
    return value


#rs scaling radius
#c concentration parameter
#h hubble constant

G = 4.301*10**(-9) #km^2 Mpc M_sun^-1 s^-2
h0 = 67 #km s-1 Mpc-1
m_tot = 2e14 # m_sun
c = 10.
omg_m = 0.3

rho_crt=3*h0**2/(8*np.pi*G)
r_200 = (3*m_tot/(4*np.pi*200*rho_crt*omg_m))**(1./3.)
rho_0 = c**3 * m_tot / (4*np.pi*r_200**3 * (np.log(1+c) - c/(1+c)))
r_s = r_200/c

print r_200
rbin = np.logspace(-2,np.log10(r_200),100)
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
mass_in = (2*np.pi*quad(lambda x: x, 0,np.min(rbin))[0])*des_pro(np.min(rbin))

print mass_in
c=0
for i in rbin:
    del_des[c] = (mass_in + 2*np.pi*quad(lambda x: x*des_pro(x), np.min(rbin), i)[0])/(np.pi*i**2) - des_pro(i)
    c=c+1




plt.xlabel('R (Mpc)')
plt.ylabel(r'$\Delta\Sigma (R)$')
plt.plot(rbin,del_des/(1e12),'r',label = 'numerical')
#plt.plot(rbin,delta(rbin),'b', label= 'interpolated')
#plt.plot(rbin,sig,'g',label = 'analytical')
plt.xlim(0.1,np.max(rbin))
#plt.ylim(0.01,101)
plt.xscale('log')
plt.yscale('log')
plt.legend()
#plt.savefig("/Users/divyarana/Desktop/model.pdf",format='pdf')
plt.show()
#'''baloo'''
