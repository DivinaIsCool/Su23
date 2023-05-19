#!/usr/bin/env python
# coding: utf-8

# ## How to import

# In[1]:


from colossus.cosmology import cosmology


# ## Set a cosmology that gives values to parameters like Omega values
# Note that length is given in Mpc/h, Mass is Msun, k is comoving h/Mpc, Time is Gyr, and Density is Msun*h^2*kpc^-3

# In[2]:


cosmo = cosmology.setCosmology('planck18')
print(cosmo)


# ## Define your own cosmology

# In[3]:


my_cosmo = {'flat': True, 'H0': 72.0, 'Om0': 0.25, 'Ob0': 0.043, 'sigma8': 0.8, 'ns': 0.97}
# use cosmo = cosmology.setCosmology('my_cosmo', **my_cosmo) to set


# ## Plot densities as function of redshift

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[5]:


zp1 = 10**np.arange(-1.0, 2.0, 0.05)
z = zp1 - 1.0

O_b = cosmo.Ob(z)
O_dm = cosmo.Om(z) - O_b
O_de = cosmo.Ode(z)
O_gamma = cosmo.Ogamma(z)
O_nu = cosmo.Onu(z)

plt.figure()
plt.loglog()
plt.xlabel('z+1')
plt.ylabel('Fraction of total density')
plt.plot(zp1, O_dm, '-', label = 'Dark matter')
plt.plot(zp1, O_b, '-', label = 'Baryons')
plt.plot(zp1, O_de, '-', label = 'Dark Energy')
plt.plot(zp1, O_gamma, '-', label = 'Photos')
plt.plot(zp1, O_nu, '-', label = 'Neutrinos')
plt.legend(loc = 4);


# ## Take inverses and derivatives

# In[6]:


z = np.array([0, 1, 2, 3])
t = cosmo.age(z)
cosmo.age(t, inverse = True, derivative = 1)


# ## Power spectrum

# In[7]:


k = 10**np.linspace(-5.0, 2.0, 500)
Pk = cosmo.matterPowerSpectrum(k)

plt.figure()
plt.loglog()
plt.xlabel('k(h / Mpc)')
plt.ylabel('P(k)')
plt.plot(k, Pk, '-');

#gives matter distribution amplitude as a function of scale


# Or its derivative

# In[8]:


Pk_deriv = cosmo.matterPowerSpectrum(k, derivative = True)

plt.figure()
plt.xscale('log')
plt.xlabel('k(h/Mpc)')
plt.ylabel('d log(P) / d log(k)')
plt.plot(k, Pk_deriv, '-');


# ## Integral quantities: correlation and variance

# In[9]:


# Correlation gives information about distribution of mass as a function of scale
R = 10**np.linspace(0.0, 2.4, 500)
xi = cosmo.correlationFunction(R, 0.0)

plt.figure()
plt.loglog()
plt.xlabel('R(Mpc/h)')
plt.ylabel('abs(xi)')
plt.plot(R, np.abs(xi), '-');
#more info can be fouend here https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html under correlation function


# In[10]:


# Variance function gives the variance of the filtered density field as a function of Radius (how much it differs from the mean)
R = 10**np.arange(0,2.4,0.005)
sigma_tophat = cosmo.sigma(R, 0.0)
sigma_sharpk = cosmo.sigma(R, 0.0, filt = 'sharp-k')
sigma_gaussian = cosmo.sigma(R, 0.0, filt = 'gaussian')

plt.figure()
plt.loglog()
plt.xlabel('R(Mpc/h)')
plt.ylabel('sigma(R)')
plt.plot(R, sigma_tophat, '-', label = 'tophat')
plt.plot(R, sigma_sharpk, '--', label = 'sharp-k')
plt.plot(R, sigma_gaussian, ':', label = 'gaussian')
plt.legend();
print(plt.ylim())
#sigma function computes the variance given by:
    # sigma^2(R, z) = (1/2(pi^2))*quad((k^2*P(k, z)*abs(W(k*R)^2)*dk), 0, infinity)
    # W is the filter functionn that can be either tophat, sharpk, or gaussian
    # P is the power spectrum defined for an array of k values and a set z value
    # sigma must be passed an array of R values and a z value


# ## Mass function
# Quantifies how many halos there are of a given mass.

# In[11]:


from colossus.lss import mass_function

z = [0.0, 1.0, 2.0, 4.0]
M = 10**np.arange(11.0, 15.5, 0.1)

plt.figure()
plt.xlabel('M200m')
plt.ylabel('dn/dln(M)')
plt.loglog()
plt.xlim(1E11, 4E15)
plt.ylim(1E-7, 1E-1)
for i in range(len(z)):
    mfunc = mass_function.massFunction(M, z[i], mdef = '200m', model = 'tinker08', q_out = 'dndlnM')
    plt.plot(M, mfunc, '-', label = 'z = %.1f' % (z[i]))
plt.legend();
#Halo mass function is given by:
  # (dn / dln(M)) = f(sigma) * (row0 / M) * |(dln(sigma**-1) / dln(M))|
  # M is a numpy array of masses, sigma is the variance on the lagrangian size scale of the halo mass in question, and f(sigma)
  # is mass function, but q_out argument gives that we want the output in dn/dlnM instead 
    
#Question: what is meant by "lagrangian size scale of the halo mass in question"?
    
# the 200m definition is one where the mass is 200 x the mean matter density


# In[12]:


z = 0.0
M = 10**np.arange(9.0, 15.5, 0.1)

plt.figure()
plt.xscale('log')
plt.xlabel('$M_{FOF}$')
plt.ylabel('$f/f_{sheth99} - 1$')
plt.xlim(1E5, 5E15)
plt.ylim(-0.6, 0.6)
plt.axhline(1.0, ls = '--', color = 'gray')

ref = mass_function.massFunction(M, z, model = 'sheth99')
for model in mass_function.models:
    if 'fof' in mass_function.models[model].mdefs:
        mfunc = mass_function.massFunction(M, z, mdef = 'fof', model = model)
        plt.plot(M, mfunc / ref - 1.0, '-', label = model)
plt.legend(loc = 3);
#this code generates a plot showing the fractional difference in the abundance of
#halos compared to a reference mass function, for different models, as a function of the halo masses.
#mass definitions include <int>m, <int>c, vir (for overdensity that varies with redshift), * (arbitrary spherical overdensity definition), fof
#Question: What is fof (friend of friend) mass definition?


# note: seppi20 model can compute conditional mass function with additional dependencies on the offset parameter (displacement between halo center of mass and the peak of its mass profile) and the Peebles spin parameter. It can also integrate over the mass (or variance) and produce the conditional mass function in spin-offset space.
# 
# Question: do we use seppi20 model?

# ## Spherical overdensity mass and radius

# In[13]:


from colossus.halo import mass_so


# In[14]:


z = np.linspace(0.0, 3.0, 40)

plt.figure()
plt.yscale('log')
plt.xlabel('z')
plt.ylabel('density (Msun h2 / kpc3)')
plt.plot(z, mass_so.densityThreshold(z, 'vir'), label = 'vir');
plt.plot(z, mass_so.densityThreshold(z, '180m'), label = '180m');
plt.plot(z, mass_so.densityThreshold(z, '180c'), label = '180c');
plt.legend();
# At high redshift, the virial overdensity threshold is equal to 180  times the matter density, and the matter densities are the same.
# At low redshift, they diverge
# Question: Why are the matter and critical densities equal at high redshift? 


# In[15]:


R = mass_so.M_to_R(1E12, 0.5, 'vir')
print(R)
M2 = mass_so.R_to_M(R, 0.5, 'vir')
print(M2)
#Get the virial radius of a halo of Mvir 10^12Msun/h at z=.5 and convert it back


# In[16]:


from colossus.halo import mass_defs
from colossus.halo import mass_adv

Mvir = 1E12
cvir = 7.0
M200c, R200c, c200c = mass_defs.changeMassDefinition(Mvir, cvir, 0.5, 'vir', '200c')
print(M200c / Mvir)
print(c200c / cvir)
#Convert from virial mass definition to 200c mass definition.
#Function requires (Mass, c, redshift, convert from definition, convert to definition)
M200c, R200c, c200c = mass_adv.changeMassDefinitionCModel(Mvir, 0.5, 'vir', '200c')
print(M200c / Mvir)
print(c200c / cvir)
#can also be used to convert by using a default model to get the concentration. remember c = Rvir / r_-2 where r_-2 is
# the radius where the log slope of the density profile is -2


# In[17]:


Mvir_z1 = 1E12
cvir_z1 = 7.0
Mvir_z0, Rvir_z0, cvir_z0 = mass_defs.pseudoEvolve(Mvir_z1, cvir_z1, 'vir', 1.0, 0.0) 
print(Mvir_z0 / Mvir_z1)
# Can see how Mvir pseudo evolves due to changing redshift. Higher redshift corresponds to lower Mvir since gravitational
# evolution creates greater overdensities over time.


# ## Halo density profiles

# In[18]:


#plot nfw density profile at z = 0 as r vs rho_nfw / rho_m
from colossus.halo import profile_nfw

Mvir = 1E15
cvir = 5.0
z = 0.0
p_nfw = profile_nfw.NFWProfile(M = Mvir, c = cvir, z = z, mdef = "vir")

r = 10**np.arange(0, 4, .02)
rho_m = cosmo.rho_m(z)
rho_nfw = p_nfw.density(r)

plt.figure()
plt.loglog()
plt.xlabel('r(kpc/h)')
plt.ylabel('density / mean')
plt.plot(r, rho_nfw / rho_m, '-', label = 'NFW');
plt.ylim(1E0, 1E7)
plt.legend();


# In[19]:


from colossus.halo import profile_einasto
from colossus.halo import profile_diemer22

# other density profile models can be used
# given a profile, you call its density with .density, mass with .MDelta or radius with .RDelta

p_einasto = profile_einasto.EinastoProfile(M = Mvir, c = cvir, z = z, mdef = 'vir')
p_d22 = profile_diemer22.ModelAProfile(M = Mvir, c = cvir, z = z, mdef = 'vir')

rho_einasto = p_einasto.density(r)
rho_d22 = p_d22.density(r)

Rvir_nfw = p_nfw.RDelta(z, 'vir')
print(p_nfw.par)
rs = p_nfw.par['rs']

print(p_einasto.par) #note that some modles give more than just the density and r value (alpha value in this case)

#same plot as above but for other models

plt.figure()
plt.loglog()
plt.xlabel('r(kpc/h)')
plt.ylabel('density / mean')
plt.plot(r, rho_nfw / rho_m, '-', label = 'NFW');
plt.plot(r, rho_einasto / rho_m, '-', label = 'Einasto');
plt.plot(r, rho_d22 / rho_m, '-', label = 'D22');
plt.axvline(Rvir_nfw, ls = '--', label = 'Rvir (NFW)');
plt.axvline(rs, ls = ':', label = 'rs (NFW)');
plt.ylim(1E0, 1E7)
plt.legend();


# In[20]:


#surface density can also be called

# Question: is the surface density just the mass per square meter at the virial radius? What exactly defines the surface?
# might be good practice to try to compute this without using the surface density colossus function.

Sigma_nfw = p_nfw.surfaceDensity(r)
Sigma_einasto = p_einasto.surfaceDensity(r)

plt.figure()
plt.loglog()
plt.xlabel('r(kpc/h)')
plt.ylabel('surface density')
plt.plot(r, Sigma_nfw, '-', label = 'NFW');
plt.plot(r, Sigma_einasto, '-', label = 'Einasto');
plt.ylim(1E6, 1E10)
plt.legend();


# In[21]:


# non linear infall and contribution from other halos causes complications at large radii (outer terms must be added)

# method 1:

from colossus.halo import profile_outer

outer_term_mean = profile_outer.OuterTermMeanDensity(z = z)
p_nfw = profile_nfw.NFWProfile(M = Mvir, c = cvir, z = z, mdef = 'vir', outer_terms = [outer_term_mean])

# method 2:

from colossus.halo import profile_composite

p_nfw = profile_composite.compositeProfile('nfw', outer_names = ['mean'], M = Mvir, c = cvir, z = z, mdef = 'vir')

# we can also create a more intricate power law in-falling profile

p_inf = profile_composite.compositeProfile('nfw', outer_names = ['mean'], M = Mvir, c = cvir, z = z, mdef = 'vir',
                                              pl_delta_1 = 10.0, pl_s = 1.4)

# or use the matter-matter correlation function times a bias

# Question: What is a bias?

from colossus.lss import bias
b = bias.haloBias(Mvir, z, 'vir')
print(b)
p_bias = profile_composite.compositeProfile('nfw', outer_names = ['cf'], M = Mvir, c = cvir, z = z, mdef = 'vir', bias = b)

# plot these and their derivatives

kwargs = dict(M = Mvir, c = cvir, z = z, mdef = 'vir')
p_mean = profile_composite.compositeProfile('einasto', outer_names = ['mean'], **kwargs)
p_cf = profile_composite.compositeProfile('einasto', outer_names = ['mean', 'cf'], bias = b, **kwargs)
p_inf = profile_composite.compositeProfile('einasto', outer_names = ['mean', 'infalling'], pl_delta_1 = 10.0, pl_s = 1.4, **kwargs)

r = 10**np.arange(0,5,0.02)

plt.figure()
plt.loglog()
plt.xlabel('r(kpc/h)')
plt.ylabel('density / mean')
plt.plot(r, p_mean.density(r) / rho_m, '-', label = 'Mean density');
plt.plot(r, p_cf.density(r) / rho_m, '-', label = 'Mean+corrfunc');
plt.plot(r, p_inf.density(r) / rho_m, '-', label = 'Mean+infalling');
plt.ylim(1E0, 1E7)
plt.legend();

plt.figure()
plt.xscale('log')
plt.xlabel('r(kpc/h)')
plt.ylabel('Logarithmic slope')
plt.plot(r, p_mean.densityDerivativeLog(r), '-', label = 'Mean density');
plt.plot(r, p_cf.densityDerivativeLog(r), '-', label = 'Mean+corrfunc');
plt.plot(r, p_inf.densityDerivativeLog(r), '-', label = 'Mean+infalling');
plt.ylim(-3.5, 0.3)
plt.legend();


# ## Profile Fitting

# In[22]:


#Create some data that is off by a little bit

from colossus.halo import mass_so

Mvir = 1E12
cvir = 7.0
Rvir = mass_so.M_to_R(Mvir, z, 'vir')
p = profile_composite.compositeProfile('nfw', outer_names = ['infalling'], M = Mvir, c = cvir, z = 0.0, mdef = 'vir', 
                            pl_delta_1 = 4.0, pl_s = 1.4, pivot = 'fixed', pivot_factor = Rvir)

# Generate random scatter around the true surface density profile
r = 10**np.arange(0.1, 3.6, 0.1)
sigma_true = p.surfaceDensity(r)
np.random.seed(155)
sigma_err = np.abs(np.random.normal(0.2, 0.1, (len(r)))) * sigma_true
sigma = sigma_true.copy()
for i in range(len(r)):
    sigma[i] += np.random.normal(0.0, sigma_err[i])

plt.figure()
plt.loglog()
plt.xlim(1E0, 5E3)
plt.ylim(5E4, 1E9)
plt.xlabel('r(kpc/h)')
plt.ylabel('Surface density')
plt.plot(r, sigma_true, '-', color = 'deepskyblue', label = 'True')
plt.legend()
plt.errorbar(r, sigma, yerr = sigma_err, fmt = '.', marker = 'o', ms = 4.0, color = 'darkblue', label = 'Data')
plt.legend();


# In[23]:


# show the parameters of this profile
print(p.par)


# In[24]:


#rho_s, r_s have to do with the standard NFW profile, while the other parameters correspond to the infalling prifle. Vary
# pl_delta_1 (the normalization at Rvir), and pl_s (the slope), but not pl_zet or pl_delta_max
mask = np.array([True, True, True, True, False, False])
x_true = p.getParameterArray(mask)
print(x_true)
ini_guess = x_true * 1.2
print(ini_guess)
p.setParameterArray(ini_guess, mask = mask)
sigma_ini = p.surfaceDensity(r)

# make a bad initial guess, then use the fit function to find a fit for the parameters

dic = p.fit(r, sigma, 'Sigma', q_err = sigma_err, mask = mask)
sigma_fit = dic['q_fit']


# In[25]:


# Question: what is chi2/Ndof and why does it matter that it's close to unity?

#plot the fit and the real function

plt.figure()
plt.loglog()
plt.xlim(1E0, 5E3)
plt.ylim(8E3, 1E9)
plt.xlabel('r(kpc/h)')
plt.ylabel('Surface density')
plt.plot(r, sigma_true, '-', color = 'deepskyblue', label = 'True')
plt.plot(r, sigma_ini, ':', color = 'purple', label = 'Initial guess')
plt.plot(r, sigma_fit, '--', color = 'firebrick', lw = 1.5, label = 'Fit')
plt.legend()
plt.errorbar(r, sigma, yerr = sigma_err, fmt = '.', marker = 'o', ms = 4.0, color = 'darkblue', label = 'Data')
plt.legend();

#note: this used a least squares fit, but an alternative option called MCMC filter can also be used

#Question: do we just use least squares, or is MCMC used?


# In[26]:


# test out a spline profile

#Question: what is a spline profile and why is it useful?

from colossus.halo import mass_so
from colossus.halo import profile_spline

#Create an NFW profile
mdef = 'vir'
p_nfw = profile_nfw.NFWProfile(M = 1E12, c = 10.0, mdef = mdef, z = 0.0)

# Create a coarse array of radii and evaluate the NFW density
Rvir = mass_so.M_to_R(Mvir, z, mdef)
r = 10**np.arange(-2.0, 1.0, 0.5) * Rvir
rho = p_nfw.density(r)

# Create a spline profile from the radius array
print('Initializing spline profile from %d radial bins.' % (len(r)))
p_spline = profile_spline.SplineProfile(r, rho = rho)

# Now create a much finer array of radii and test how well the spline does
r = 10**np.arange(-2.0, 1.0, 0.01) * Rvir
rho_m = cosmo.rho_m(z)
rho_nfw = p_nfw.density(r)
rho_spline = p_spline.density(r)

plt.figure()
plt.xscale('log')
plt.xlabel('r(kpc/h)')
plt.ylabel('spline / true')
plt.plot(r, rho_spline / rho_nfw, '-');

#Notice that spline profiles are not as useful at high radii

#Also note that we can create our own density profiles using classes (shown at the end of the profiles tutorial)

#Question: Do we make our own profiles, or use existings ones? If we use existing ones, which ones do we typically use?


# ## Concentration
# halo concentration is defined as the ratio fo the outer radius of a halo to the scale radius

# In[27]:


from colossus.halo import concentration

# get list of models

for model_name in concentration.models:
    print(model_name)


# In[28]:


# find concentrations of given masses at given redshift with a certain model

M = 10**np.arange(10.0, 15.0, 1.0)
concentration.concentration(M, 'vir', 0.5, model = 'bullock01')

#some masses or redshifts can be outside the validity of the concentration model, use the following method to verify validity:

c, mask = concentration.concentration(M, 'vir', 0.5, model = 'dutton14', range_return = True)
print(mask)
print(c[mask])


# ## Splashback radius
# splashback radius is where the halo density profile undergoes a significant drop, where the logarithmic slope experiences a sudden change.

# In[29]:


#list the varius splashback models
#note that different authors mean different things by "splashback radius"
#Question: do we use the splashback radius, and if so what model / models? How do we typically define it?
from colossus.halo import splashback

for model_name in splashback.models:
    print(model_name)


# In[30]:


#plot of splashback radii against accretion rate for different models
#note that the peak density was needed
#Question: Why exactly is the splashback radius associated with the accretion rate of the halo?

from colossus.lss import peaks

z = 0.0
M200m = 1E12
nu200m = peaks.peakHeight(M200m, z)
Gamma = np.arange(0.0, 5.1, 0.1)

plt.figure()
plt.xlabel('Accretion rate')
plt.ylabel('Rsp/R200m')
for model_name in splashback.models:
    RspR200m, mask = splashback.splashbackModel('RspR200m', Gamma = Gamma, nu200m = nu200m, z = z, 
                                model = model_name, rspdef = 'sp-apr-p75', statistic = 'median')
    plt.plot(Gamma[mask], RspR200m, label = model_name.replace('_', '\_'))
plt.legend();


# In[31]:


# the splashbackRadius function can be used to simply calculated Rsp or Msp
# this can be done without a mass accretion rate, but including one will make the measurement more accurate

z = 0.0
mdef = 'vir'
Mvir = 1E12
cvir = 10.0

Rsp, Msp, _ = splashback.splashbackRadius(z, mdef, M = Mvir, c = cvir, rspdef = 'sp-apr-p90')
print(Rsp, Msp)

Rsp, Msp, _ = splashback.splashbackRadius(z, mdef, M = Mvir, c = cvir, Gamma = 3.0, rspdef = 'sp-apr-p90')
print(Rsp, Msp)


# In[ ]:




