#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import statement:
#note that there was a glitch with matplotlib and I had to use that command to reset the matplotlib settings

import disSat as dis

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update(plt.rcParamsDefault)
get_ipython().run_line_magic('matplotlib', 'inline')

from colossus.cosmology import cosmology
from colossus.halo import mass_so
from colossus.halo import mass_defs
from colossus.halo import mass_adv
from colossus.halo import profile_nfw
from colossus.halo import profile_einasto
from colossus.halo import concentration
from scipy.integrate import quad
from scipy.integrate import quadrature

cosmo = cosmology.setCosmology('planck18')


# In[49]:


# Generate Mass array and compute R, c, and r_s for later use

scatter_smhm = 0.15
scatter_concentration = 0.16
scatter_rhalf2D = 0.234

n = 100000
M200 = 10**np.linspace(7, 12, n)
R200 = mass_so.M_to_R(M200*cosmo.h, 1, '200c') / cosmo.h
c200 = np.random.lognormal(mean=np.log(concentration.concentration(M200*cosmo.h, '200c', 1, 'diemer19')), 
                           sigma=(scatter_concentration / np.log(10)))
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

m_stellar = np.random.lognormal(mean=np.log((M200*2*N)*((((M200/M1)**(-1))*Beta)+((M200/M1)**Gamma))**(-1)), sigma=(scatter_smhm / np.log(10)))
Reff = np.random.lognormal(mean=np.log(10**(.268*np.log10(m_stellar)-2.11)), sigma=(scatter_rhalf2D / np.log(10)))
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

nbins = 5000
sigma_split = np.array_split(sigma_los_NFW, nbins)
quant_max = np.array([])
quant_min = np.array([])
for sub_array in sigma_split:
    qmax = np.quantile(sub_array, q=(0.5 + (0.68/2)))
    qmin = np.quantile(sub_array, q=(0.5 - (0.68/2)))
    quant_max = np.append(quant_max, qmax)
    quant_min = np.append(quant_min, qmin)
downsampled_Reff = Reff[::20] 

# Plot sigma_los vs R_eff (with scatter)

plt.figure()
plt.loglog()
plt.ylabel('(<sigma_los^2)^1/2>(km/s)')
plt.xlabel('R_eff')
#plt.plot(Reff * 1000, sigma_los_NFW, '-', label = 'NFW');
plt.fill_between(downsampled_Reff * 1000, quant_min, quant_max, alpha=0.3, label='NFW')
#plt.plot(M200, sigma_los_CORED, '-', label = 'cored');
plt.gca().set_yticks([1e1])
plt.gca().set_xticks([1e2, 1e3])
plt.xlim(15, 4000)
plt.ylim(1.5, 23)
plt.legend();


# In[ ]:




