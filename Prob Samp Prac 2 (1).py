#!/usr/bin/env python
# coding: utf-8

# The assignment I was given:
# 
# You could move forward with probabilistic sampling. Assuming that:
# 1. The probability for a halo to have mass M is proportional to M^(-2}
# 2. The probability of a halo with mass M to have stellar mass M_star follows a lognormal distribution with mean=0.05 * M * (M/1e12 M_sun) and sigma=0.35
# 3. The probability for a halo with stellar mass M_star to have a half-light radius R_eff follows a lognormal distribution with mean=Eq. (5) in arXiv:2106.09050 (with the correction I mentioned a few days ago) and sigma=0.54
# 
# Could you generate a 2D histogram of R_eff and M_star of halos with masses between 1e7 M_sun and 1e10 M_sun? Functions you may find useful:
# 1. np.random.choice
# 2. np.random.lognormal
# 3. plt.hist2d

# In[2]:


#import numpy and matplotlib

import numpy as np
import matplotlib.pyplot as plt

#define initial constants like sample numberes and max and min masses to be evaluated over

nsamp = 10000
min_mass_halo = 1e7
max_mass_halo = 1e10

#define the halo masses and sample them to obtain an array of sampled halo masses

halo_masses_tobesampled = np.linspace(min_mass_halo, max_mass_halo, nsamp)
prob_halo = halo_masses_tobesampled**(-2)
prob_halo /= np.sum(prob_halo) #gotta normalize because the probabilities have to add up to 1
halo_samples = np.random.choice(halo_masses_tobesampled, nsamp, p=prob_halo)

#compute the array of stellar masses

SM_sigma = 0.35
SM_mean = .05 * halo_samples * (halo_samples / 1e12)
M_star = np.random.lognormal(mean=(SM_mean), sigma=SM_sigma)

#compute the array of sampled R_eff

Reff_sigma = 0.54
Reff_mean = 10**(0.268 * np.log10(M_star) - 2.11)
R_eff = np.random.lognormal(mean=(Reff_mean), sigma=Reff_sigma)

#plot the 2d histogram of stellar masses vs R_eff

plt.hist2d(R_eff, M_star, bins=230, range=[[0, 6], [0, 5000]], cmap='hot')
plt.colorbar(label='Count')
plt.xlabel('2D Half-Light Radius (kpc)')
plt.ylabel('Stellar Mass (Msun)')
plt.title('2D Histogram of R_eff and M_star for Halos')
plt.show()
plt.hist2d(R_eff, M_star, bins=200, range=[[0, 3], [0, 800]], cmap='hot')
plt.colorbar(label='Count')
plt.xlabel('2D Half-Light Radius (kpc)')
plt.ylabel('Stellar Mass (Msun)')
plt.title('Zoomed In')
plt.show()


# In[3]:


print(np.random.lognormal(mean=(SM_mean), sigma=SM_sigma))


# In[ ]:




