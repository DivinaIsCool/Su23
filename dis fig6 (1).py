#!/usr/bin/env python
# coding: utf-8

# In[83]:


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

from disSat.vutils import menc
from disSat.vutils import mvir2sigLOS
from disSat.vutils import Reff as R_eff

cosmo = cosmology.setCosmology('planck18')


# In[84]:


# LEFT SIDE NO SCATTER

n=2500
M200 = 10**np.linspace(7, 12, n)
Reff_NFW = R_eff(m200=M200, profile='NFW', zin=1)
sigLOS_NFW = mvir2sigLOS(mvir=M200, profile='NFW')
Reff_coreNFW = R_eff(m200=M200, profile='coreNFW', zin=1)
sigLOS_coreNFW = mvir2sigLOS(mvir=M200, profile='coreNFW')
Reff_SIDM = R_eff(m200=M200, profile='SIDM', sigmaSI = 1, zin=1)
sigLOS_SIDM = mvir2sigLOS(mvir=M200, profile='SIDM', sigmaSI = 1)
plt.figure()
plt.loglog()
plt.xlabel('Reff')
plt.ylabel('(<sigma_los^2)^1/2>(km/s)')
plt.plot(Reff_NFW*1000, sigLOS_NFW, '-', label = 'NFW');
plt.plot(Reff_coreNFW*1000, sigLOS_coreNFW, '-', label = 'coreNFW');
plt.plot(Reff_SIDM*1000, sigLOS_SIDM, '-', label = 'SIDM');
plt.gca().set_yticks([1e1])
plt.gca().set_xticks([1e2, 1e3])
plt.xlim(15, 4000)
plt.ylim(1.5, 23)
plt.legend();


# In[111]:


# LEFT SIDE WITH SCATTER

#define the scatter constants

scatter_smhm = 0.15
scatter_concentration = 0.16
scatter_rhalf2D = 0.234

#define mass range

n=100000
M200 = 10**np.linspace(7, 12, n)
#M200 = np.logspace(7, 12, n)
#compute concentration

c200 = np.random.lognormal(mean=np.log(cNFW(M200, z=1, massdef='200c')), sigma=(scatter_concentration / np.log(10)))

#compute stellar mass

m_stellar = np.random.lognormal(mean=np.log(moster13(mhalo=M200, z=1)), sigma=(scatter_smhm / np.log(10)))

#compute Reff and sigLOS for each model

Reff_NFW = np.random.lognormal(mean=np.log(R_eff(M200, profile='NFW', mstar = m_stellar, c200=c200, zin=1)), sigma=(scatter_rhalf2D / np.log(10)))
sigLOS_NFW = mvir2sigLOS(M200, profile='NFW', c200=c200, mstar=m_stellar)

Reff_coreNFW = np.random.lognormal(mean=np.log(R_eff(M200, profile='coreNFW', mstar = m_stellar, c200=c200, zin=1)), sigma=(scatter_rhalf2D / np.log(10)))
sigLOS_coreNFW = mvir2sigLOS(M200, profile='coreNFW', c200=c200, mstar=m_stellar)

Reff_SIDM = np.random.lognormal(mean=np.log(R_eff(M200, profile='SIDM', mstar = m_stellar, c200=c200, zin=1, sigmaSI=1)), sigma=(scatter_rhalf2D / np.log(10)))
sigLOS_SIDM = mvir2sigLOS(M200, profile='SIDM', c200=c200, mstar=m_stellar, sigmaSI=1)

#plot sigLOS vs Reff for each model

plt.figure()
plt.loglog()
plt.xlabel('Reff (pc)')
plt.ylabel('(<sigma_los^2)^1/2> (km/s)')

plt.plot(Reff_NFW*1000, sigLOS_NFW, '-', label = 'NFW', alpha=.5);

plt.plot(Reff_coreNFW*1000, sigLOS_coreNFW, '-', label = 'coreNFW', alpha=.5);

plt.plot(Reff_SIDM*1000, sigLOS_SIDM, '-', label = 'SIDM', alpha=.5);

plt.gca().set_yticks([1e1])
plt.gca().set_xticks([1e2, 1e3])
plt.xlim(15, 4000)
plt.ylim(1.5, 23)

plt.legend();


# In[110]:


# RIGHT SIDE 

n = 100000
M200 = 10**np.linspace(6, 12, n)

c200 = np.random.lognormal(mean=np.log(cNFW(M200, z=1, massdef='200c')), sigma=(scatter_concentration / np.log(10)))

m_stellar = np.random.lognormal(mean=np.log(moster13(mhalo=M200, z=1)), sigma=(scatter_smhm / np.log(10)))

Reff_NFW = np.random.lognormal(mean=np.log(R_eff(M200, profile='NFW', mstar = m_stellar, c200=c200, zin=1)), sigma=(scatter_rhalf2D / np.log(10)))
Meff_NFW = menc(Reff_NFW, M200, profile='NFW', zin=1, mstar=m_stellar, c200=c200)

Reff_coreNFW = np.random.lognormal(mean=np.log(R_eff(M200, profile='coreNFW', mstar = m_stellar, c200=c200, zin=1)), sigma=(scatter_rhalf2D / np.log(10)))
Meff_coreNFW = menc(Reff_coreNFW, M200, profile='coreNFW', zin=1, mstar=m_stellar, c200=c200)

Reff_SIDM = np.random.lognormal(mean=np.log(R_eff(M200, profile='SIDM', mstar = m_stellar, c200=c200, zin=1, sigmaSI=1)), sigma=(scatter_rhalf2D / np.log(10)))
Meff_SIDM = menc(Reff_SIDM, M200, profile='SIDM', zin=1, mstar=m_stellar, c200=c200, sigmaSI=1)

L = m_stellar / 2
M = -2.5*np.log10(L)

mL_NFW = Meff_NFW / L
mL_coreNFW = Meff_coreNFW / L
mL_SIDM = Meff_SIDM / L

plt.figure()
plt.yscale('log')
plt.xlabel('M_v')
plt.ylabel('Dynamical mass-to-light ratio < Reff (Msun/Lsun)')

plt.plot(M, mL_NFW, '-', label = 'NFW', alpha=.5);

plt.plot(M, mL_coreNFW, '-', label = 'coreNFW', alpha=.5);

plt.plot(M, mL_SIDM, '-', label = 'SIDM', alpha=.5);

plt.gca().set_yticks([1e1, 1e2, 1e3, 1e4])
plt.gca().set_xticks([-2.5, -5, -7.5, -10, -12.5, -15])
plt.xlim(0, -17.5)
plt.ylim(1, 1e4)

plt.legend();


# In[ ]:




