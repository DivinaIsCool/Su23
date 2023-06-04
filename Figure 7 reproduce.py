#!/usr/bin/env python
# coding: utf-8

# In[4]:


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


# In[74]:


npops = 100
mw_satpops_wdm_2 = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_wdm_2:
    satpop.dark_matter = dis.dark_matter.models.WDM(mWDM=2.)
    satpop.generate_population()
    
mw_satpops_wdm_4 = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_wdm_4:
    satpop.dark_matter = dis.dark_matter.models.WDM(mWDM=4.)
    satpop.generate_population()

mw_satpops_wdm_6 = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_wdm_6:
    satpop.dark_matter = dis.dark_matter.models.WDM(mWDM=6.)
    satpop.generate_population()    
    
mw_satpops_wdm_12 = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_wdm_12:
    satpop.dark_matter = dis.dark_matter.models.WDM(mWDM=12.)
    satpop.generate_population()    

mw_satpops_CDM = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_CDM: 
    satpop.density_profile = 'nfw'
    satpop.generate_population()

dis.plot.plot_theoretical_velocity_fuction(mw_satpops_CDM, label='CDM, NFW', color='black')
dis.plot.plot_theoretical_velocity_fuction(mw_satpops_wdm_12, label='WDM, 12 keV', color='C0')    
dis.plot.plot_theoretical_velocity_fuction(mw_satpops_wdm_6, label='WDM, 6 keV', color='C1')
dis.plot.plot_theoretical_velocity_fuction(mw_satpops_wdm_4, label='WDM, 4 keV', color='C2')
dis.plot.plot_theoretical_velocity_fuction(mw_satpops_wdm_2, label='WDM, 2 keV', color='C3')
dis.plot.plot_observed_velocity_function()
dis.plot.plot_corrected_velocity_function(verbose=False)
dis.plot.finalize(legend=True)


# In[22]:


npops = 100
mw_satpops_wdm_2 = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_wdm_2:
    satpop.dark_matter = dis.dark_matter.models.WDM(mWDM=2.)
    satpop.generate_population()
    
dis.plot.plot_theoretical_velocity_fuction(mw_satpops_wdm_2, label='WDM, 2 keV', color='C1', alpha = 1)
dis.plot.finalize(legend=True)


# In[ ]:


mw_satpops_CDM = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_CDM: 
    satpop.density_profile = 'nfw'
    satpop.generate_population()
mw_satpops_CDM_ml1 = [ dis.satpops.MilkyWaySatellites() for i in range(npops) ]
for satpop in mw_satpops_CDM_ml1: 
    satpop.density_profile = 'nfw'
    satpop.mleft = 0.1
    satpop.generate_population()
dis.plot.plot_theoretical_velocity_fuction(mw_satpops_CDM_ml1, label=False, color='black', linestyle='.')
dis.plot.plot_theoretical_velocity_fuction(mw_satpops_CDM, label='CDM', color='black')
dis.plot.finalize(legend=True)


# In[54]:


dis.plot.plot_corrected_velocity_function(verbose=False)
dis.plot.finalize(legend=True)


# In[ ]:




