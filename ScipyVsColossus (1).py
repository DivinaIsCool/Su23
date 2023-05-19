#!/usr/bin/env python
# coding: utf-8

# ## Attempt to use scipy to compute sigma variance and compare to variance computed using colossus

# In[1]:


#import necessary modules
import matplotlib.pyplot as plt
import numpy as np
import scipy
from colossus.cosmology import cosmology
from scipy.integrate import quad


# In[2]:


#set the cosmology as planck18
cosmo = cosmology.setCosmology('planck18', persistence = 'r')


# In[21]:


#create a matter power spectrum and k range to be used by both scipy and colossus
k = 10**np.linspace(-5.0, 2.0, 500)
Pk = cosmo.matterPowerSpectrum(k)
plt.figure()
plt.loglog()
plt.xlabel('k')
plt.ylabel('Pk')
plt.plot(k, Pk, '-', label = 'tophat scipy')


# In[25]:


#create R range for the sigma function to be used in the scipy method and colossus method
R = 10**np.linspace(0,2.4,500)
#define Heaviside function for sharp-k filter (this is what Chat GPT told me it is)
def heaviside(x):
    return np.where(x < 0, 0, 1)
#define the W function three ways to be used in the scipy method
def Wk_gaus(k, R): 
    return (np.exp((-(k*R)**2)/2)) * R/np.sqrt(2*np.pi)
def Wk_tophat(k, R):
    return ((3)/(k*R)**3)*(np.sin(k*R)-k*R*np.cos(k*R))
def Wk_sharpk(k, R):
    return heaviside(1 - k*R)
#make a function for the integrand
def integrand(k, Pk, Wkfunction, R):
    return k**2 * cosmo.matterPowerSpectrum(k) * np.absolute(Wkfunction(k, R))**2
#create sigma function and compute it for each Window function
def sigma(Wfunc):
    integral_gaus_values = np.array([])
    for i in range(len(R)):
        integral_gaus_value, gauserr = quad(integrand, 0, 10000, args=(Pk, Wk_sharpk, R[i]))
        integral_gaus_values = np.append(integral_gaus_values, [integral_gaus_value])
    return np.sqrt((1/2*(np.pi)**2) * integral_gaus_values)
#sigma_gaus_scipy = sigma(Wk_gaus(k, R))
#sigma_tophat_scipy = sigma(Wk_tophat(k, R))
sigma_sharpk_scipy = sigma(Wk_sharpk(k, R))

#Question: What is heaviside function?


# In[5]:


#compute the three sigma functions with colossus
sigma_tophat_colossus = cosmo.sigma(R, 0.0)
sigma_sharpk_colossus = cosmo.sigma(R, 0.0, filt = 'sharp-k')
sigma_gaus_colossus = cosmo.sigma(R, 0.0, filt = 'gaussian')


# In[28]:


#plot everything
plt.figure()
plt.loglog()
plt.xlabel('R(Mpc/h)')
plt.ylabel('sigma(R)')
plt.plot(R, sigma_tophat_scipy, '-', label = 'tophat scipy')
plt.plot(R, sigma_sharpk_scipy, '--', label = 'sharp-k scipy')
plt.plot(R, sigma_gaus_scipy, ':', label = 'gaussian scipy')
plt.plot(R, sigma_tophat_colossus, '-', label = 'tophat colossus')
plt.plot(R, sigma_sharpk_colossus, '--', label = 'sharp-k colossus')
plt.plot(R, sigma_gaus_colossus, ':', label = 'gaussian colossus')
#plt.ylim(0.0023585901002120147, 3.3813731968470884)
plt.legend();

#Question: Why is it so off? Does it have to do with the way the integrand is defined? 
# Dodleson defines it different from the colossus documentation.


# In[ ]:




