#!/usr/bin/env python
# coding: utf-8

# In[6]:


# Goal is to reproduce NFW and cored (mleft=1) in fig 4 in https://arxiv.org/pdf/2106.09050.pdf

from colossus.cosmology import cosmology
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from colossus.halo import mass_so
from colossus.halo import mass_defs
from colossus.halo import mass_adv
from colossus.halo import profile_nfw
from colossus.halo import profile_einasto
from colossus.halo import concentration
from scipy.integrate import quad

cosmo = cosmology.setCosmology('planck18')

#They start with Mvir=M200=1e9Msun. Derive R200, c200, r_s, and generate a list of radial values for the plot

M200 = 1E9
R200 = mass_so.M_to_R(M200, 1, '200c')
c200 = concentration.concentration(M200, '200c', 1, 'diemer19')
r_s = R200/c200
r = np.linspace(1e-3, R200, 1500)

# Define density functions and integrate out to R200 to obtain enclosed mass functions

def rho_prof_NFW(r, C):
    return C / ((r/r_s)*((1+r/r_s)**2))

def rho_prof_CORED(r, C):
    return C / ((1+(r/r_s)**3))

result1, _ = quad(lambda r: 4*np.pi*rho_prof_NFW(r, 1)*r**2, 1e-3, R200)
C1 = M200/result1

result2, _ = quad(lambda r: 4*np.pi*rho_prof_CORED(r, 1)*r**2, 1e-3, R200)
C2 = M200/result2

Menc_NFW = np.array([])
for i in r:
    result1, _ = quad(lambda r: 4*np.pi*rho_prof_NFW(r, C1)*r**2, 1e-3, i)
    Menc_NFW = np.append(Menc_NFW, result1)
    
Menc_CORED = np.array([])
for i in r:
    result2, _ = quad(lambda r: 4*np.pi*rho_prof_CORED(r, C2)*r**2, 1e-3, i)
    Menc_CORED = np.append(Menc_CORED, result2)

# Plot enclosed mass vs distance from center     
    
plt.figure()
plt.loglog()
plt.xlabel('r/r200')
plt.ylabel('M(<r) / Mtot,0')
plt.plot(r / R200, Menc_NFW / M200, '-', label = 'NFW');
plt.plot(r / R200, Menc_CORED / M200, '-', label = 'cored');
plt.gca().set_yticks([1e-7, 1e-5, 1e-3, 1e-1])
plt.xlim(1e-3, R200/8)
plt.legend();


# In[5]:


# Generate Mass array and compute R, c, and r_s for later use

M200 = 10**np.linspace(8, 11, 2500)
R200 = mass_so.M_to_R(M200*cosmo.h, 1, '200c') / cosmo.h
c200 = concentration.concentration(M200*cosmo.h, '200c', 1, 'diemer19')
r_s = R200/c200

# Define constants and compute the important constants for the SHM relation

M10 = 11.59
M11 = 1.195
N10 = .0351
N11 = -.0247
Beta10 = 1.376
Beta11 = -0.826
Gamma10 = 0.608
Gamma11 = 0.329
M1 = 10**(M10 + .5*M11)
N = N10 + .5*N11
Beta = Beta10 +.5*Beta11
Gamma = Gamma10 + .5*Gamma11

# Find stellar mass from SHM relation, compute Reff from that and r_half from that

m_stellar = (M200*2*N)*(((M200/M1)**(-1*Beta))+(M200/M1)**Gamma)**(-1)
Reff = 10**(.268*np.log10(m_stellar)-2.11)
r_half = Reff / .75

# Integrate along density profile from 1e-3 to r_half for each value to obtain M_half array for each model
# Note that this requires a unique C value for each integral, so more loops were necessary to obtain those

def rho_prof_NFW(r, rs, C):
    return C / ((r/rs)*((1+r/rs)**2))

def rho_prof_CORED(r, rs, C):
    return C / ((1+(r/rs))**3)

Cvalues1 = np.array([])
for i, o in enumerate(R200):
    result1, _ = quad(lambda r: 4*np.pi*rho_prof_NFW(r, r_s[i], 1)*r**2, 1e-3, o)
    C1 = M200[i]/result1
    Cvalues1 = np.append(Cvalues1, C1)

Cvalues2 = np.array([])    
for i, o in enumerate(R200):
    result2, _ = quad(lambda r: 4*np.pi*rho_prof_CORED(r, r_s[i], 1)*r**2, 1e-3, o)
    C2 = M200[i]/result2
    Cvalues2 = np.append(Cvalues2, C2)

M_half_NFW = np.array([])
for i, o in enumerate(r_half):
    result, _ = quad(lambda r: 4*np.pi*rho_prof_NFW(r, r_s[i], Cvalues1[i])*r**2, 1e-3, o)
    M_half_NFW = np.append(M_half_NFW, result)

M_half_CORED = np.array([])
for i, o in enumerate(r_half):
    result, _ = quad(lambda r: 4*np.pi*rho_prof_CORED(r, r_s[i], Cvalues2[i])*r**2, 1e-3, o)
    M_half_CORED = np.append(M_half_CORED, result)

# Compute sigma_los    
    
G = 4.302 * 10**-6   # units of km^2*kpc*s^-2*Msun^-1
sigma_los_NFW = np.sqrt((G*(M_half_NFW))/(4*Reff))
sigma_los_CORED = np.sqrt((G*(M_half_CORED))/(4*Reff))

# Plot sigma_los vs M200

plt.figure()
plt.loglog()
plt.xlabel('M200(Msun)')
plt.ylabel('(<sigma_los^2)^1/2>(km/s)')
plt.plot(M200, sigma_los_NFW, '-', label = 'NFW');
plt.plot(M200, sigma_los_CORED, '-', label = 'cored');
plt.gca().set_yticks([1e0, 1e1])
plt.gca().set_xticks([1e8, 1e9, 1e10, 1e11])
plt.legend();


# In[3]:


x = .234*np.random.poisson(0, 25)
print(x)


# In[ ]:





# In[ ]:




