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
m_tot = 1e16 # m_sun
c = 10.
omg_m = 0.3

rho_crt=3*h0**2/(8*np.pi*G)
r_200 = (3*m_tot/(4*np.pi*200*rho_crt*omg_m))**(1./3.)
rho_0 = c**3 * m_tot / (4*np.pi*r_200**3 * (np.log(1+c) - c/(1+c)))
r_s = r_200/c

print r_200
rbin = np.logspace(-2,np.log10(r_200),30)
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

# DAUGHTER HALO PARAMETERS
# position of center of the disc
x0 = r_200/2
y0 = r_200/2
r0 = np.sqrt(x0**2 + y0**2)
rd_bin = rbin/4

del_sigma = np.zeros(len(rd_bin))

# integrating over theta getting mean over a circle
def des_cir(r,r0):
    value = quad(lambda j: des_pro(np.sqrt(r0**2 + r**2 + 2*r0*r*np.cos(j))), -np.pi, np.pi)[0]/(2*np.pi)
    return value


c = 0
for r in rd_bin:
    del_sigma[c] = (2*np.pi*quad(lambda r: des_cir(r,r0), 0,r)[0])/(np.pi*r**2) - des_cir(r,r0)
    print 'hello %2.5f' % (r)
    c = c+1
plt.xlabel('R (Mpc)')
plt.ylabel(r'$\Delta\Sigma (R)$')
plt.plot(rd_bin,del_sigma,'r',label = 'sattelite numerical')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

