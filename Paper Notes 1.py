#!/usr/bin/env python
# coding: utf-8

# ## Notes for "The Milky Way satellite velocity function is a sharp probe of small-scale
# ## structure problems" Paper
# By: Stacy Y. Kim and Annika H.G. Peter
# https://arxiv.org/pdf/2106.09050.pdf

# In[2]:


# Introduction

import numpy as np

# Two major problems: CDM predicts more satellite galaxies around the milky way by an order of magnitude, but we do not
# observe this (preior to the 2000s). Also, galaxy profiles are less dense in the center than predicted.

# After 2000s, more observations were made, and we also discovered that halos above 10**8 Msun must be populated by visible
# galaxies. Depending on one's hypothesis of the radial distribution of satellited in the MW, the milky way may have a
# "too many satellites problem"

# For the core/cusp problem, many argue that baryons play a big role in the dynamics of dark matter on far scales. Many also 
# propose that dark matter self interactions on small scales could bring observations closer to theoretical predictions. 

# The halo properties predicted do not match the observed kinematics. Origninally, to solve the missing satellites problem we
# used 

v_circ = np.sqrt(3)*sigma_los_star # where sigma_los is the line of sight stellar velocity dispersion of satellites and 
# v_circ is the halo circular velocity

# but this paper challenges this notion to find a better relationship to map the accurate amount of satellite galaxies with
# a certain v_circ

# Question: How does LOS relate to translational motion and radial motion? Is one omitted?

# Question: Why does that v_circ relation exist? 


# In[ ]:


# 2.1 Completeness Correction

from scipy.integrate import quad

# apply completeness correction:

def c(sigma_los_star, L):
    def integrand(r):
        return n(r)
    integral_top, _ = quad(integrand, V_vir)
    integral_botton, __ = quad(integrand, V_obs(L))
    return integral_top / integral_bottom

# V_obs(L) is the volume over which satellites of luminosity L have been surveyed
# V_Vir is the virial volume of the milky way
# n(r) is the 3D satellite distribution

# We estimate the total number of satellites with LOS stellar velocity dispersion sigma_los with:

def N(>sigma_los):
    sum(c(sigma_los_star, L), i=1, N_obs) #sum is a summation function with start of 1 and upper limit N_obs
    
# Question (bottom left of page 3): Why assumption that satellites follow the distribution of the MW's dark matter halo?

# When satellite galaxies are being stripped, their central regions likely retain galaxy-like morphology inside the milky way, 
# signalling that they can survive the tital forces of the MW. (assumed by NFW)

# Question: What is anisotropy? (Example, it says that a 30% scatter was added to account for anisotropy)

# Question: What is the average amount of stars used to determine a sigma_los?


# In[ ]:


# 2.2 Observational Uncertainties

# Consider sources of uncertainty in this section. Ultrafaint satellites have far fewer star measurements than normal satellites
# do. Binary stars may be driving up the velocity dispersion measurements because they appear more in those ultrafaint satellites
# Bootes II dwarf galaxy sigma_los values have a big effect on the completeness corrected VF. 
# In general, there is a lot of uncertainty with the VF because it depends so highly on the accurate measurements of sigma_los.

# Tidal stripping does not affect sigma_los until the subhalo is stripped to about 1% of its mass
# Earlier we assumed that analog galaxies with similar luminosities have similar sigma_los, but this assumption may fail if we
# take into account tital stripping. Determining which galaxies have undergone severe stripping is non-trivial.
# Galaxies do not exhibit stellar mass loss until >90% of halo mass has been stripped. Some galaxies, like Hercules exhibit this.

# Question: What does metal richness have to do with the likelihood of tidal disruption? *

# If we assume that satellites like Hercules and a few others have been stripped to 1%, the VF shifts dramatically at 10km/s


# In[ ]:


# 3 Theoretical Predictions for the Velocity Function

# This paper adopts an analytic approach to test scaling relations to get theoretical predictions for the VF

# Question: What does it mean to adopt a normalization that is lower by 20%? (Page 5 bottom left) What is normalization?

# Host halos exhibit a scatter in the number of subhalos that they have due to accretion histories. They generated 100 
# realizations of MW subhalo populations.

# Question: What does equation 3 mean? Poisson component? Intrinsic component? How does this relate to dex (page 6 top left)?

# Question: Can you elaborate on how reionization affected halo formation, and what a reionization redshift is?

# They use some conservative assumptions about v_max in pre-reionization, present day, and halo fraction that contains galaxies
# to obtain the observable satellite mass function. To compare this with the completeness-corrected VFs, they map the 
# halo masses into:

<sigma_los^2> = (G * M_(1/2)) / (4 * R_eff) # where R_eff is the 2D halflight radius (r_1/2 = (4/3)R_eff) and M_1/2 is the mass
# enclosed in the 3D halflight radius.

log_10(R_eff) = .239log_10(M_*) - 1.68 # empirical relationship between stellar mass of isolated dwarfs and R_eff

# Given a density profile for the halos, we can derive M_1/2 by integrating the profile to r_1/2. Stripping modifies halo 
# density profiles. Dwarfs with masses about equal to 10**7 Msun have baryonic components that contribute significantly to M_1/2
# To account, they computed M_* and added half of it to the dm mass within that radius (because halflight radius).


# In[ ]:


# 3.1 Density Profiles

# This paper consideres two alpha, beta, gamma broken power law profiles:
(alpha_1, beta_1, gamma_1) = (1, 3, 1)
(alpha_2, beta_2, gamma_2) = (1, 3, 0)

# Question: Does this mean that the density profile goes as rho(r) ∝ r**alpha, then ∝ r**beta, then ∝ r**gamma ?

# They also consider an SIDM model and a coreNFW model that assigns core sizes based on the length of time a dwarf has formed
# stars: t_SF. 


# In[ ]:


# 3.2 Tidal Stripping

# Tidal stripping takes mass from the galaxy and subhalo, affecting R_eff and sigma_los. 
# Tidal stripping effects on suhalos are modeled according to fitting functions derived in past simulations. This found that
# alpha,beta,gamma profiles were largely unchanged after stripping, except on outskirts, where slope becomes gamma = 5.
# coreNFW and SIDM models haven't had these studies done, so this paper uses its analytical approach to find the effects
# of tidal stripping. 
# In figure 4 and in the remainder, they show the results for unstripped halos (m_left = 1),
# 90% stripped halos (m_left = .1), and 99% stripped halos (m_left = .01).

# Stellar mass does not begin to deplete until halo mass has been stripped by >90%. Well after 99% has been stripped, stellar 
# mass will have dropped by a factor of 2. This paper considers this stripping unlikely for the vast majority of satellites.


# In[ ]:


# 3.3 Comparisons Against Observed Relations

# In figure 6 they compare their results with observations by plotting sigma_los vs R_eff and Mass/Luminosity ratio < R_eff vs 
# absolute magnitude. coreNFW best matches. SIDM with .1 cm^2/g cross section is better
# than the originally assumed 1 cm^2/g

# Why plot Mass/Luminosity ratio < R_eff vs Absolute magnitude?


# In[ ]:


# 4.1 Results: Cusps and cores

# They present the velocity functions for each profile and assuming varried levels of stripping

# Question (page 9 right middle): What does it mean to assume star formation continued until infall onto the MW at z_in = 1?

# Concerning the top left of figure 7:

# If subhalos are cuspy (1,3,1), the theoretical VF can explain the most conservative completeness-corrected
# VF as long as the average satellite’s halo is stripped by no more than 90% of its dark matter

# However, the correction assuming disk stripping requires far more satellites than the theoretical
# predictions (a too many satellites problem).

# At intermediate sigma_los ~ 10km/s we recover a missing satellites problem

# core alpha,beta,gamma profiles (1,3,0) cannot reach even the most conservative completeness-corrected VF, so this model is generally
# disfavored

# Concerning top right of figure 7:

# We find that density profiles arising from realistic models of baryonic feedback on dark-matter halos better capture the VF relative to
# the 𝛼𝛽𝛾 NFW or cored models

# We find a switch between cusps at M_200 < 5 x 10^8M_sun and cores above that value. The VF model based on these criteria
# better matches the completeness corrected VF. The exact threshold for this switch is uncertain. 


# In[ ]:


# 4.2 Results: SIDM

# Concerning bottom left of figure 7:

# SIDM cross sections of ~.5cm^2/g were more accurate at high sigma_los, but completely off at low values. Some argue that a 
# velocity dependent cross section is a valid model and can reproduce observed density patterns.


# In[ ]:


# 4.3 Results: WDM

# The satellite VF is unique probe of WDM models due to its sensitivity to the supression of power and internal densities on
# low mass scales

# Question: By "power" do they mean the power spectrum?

# As halos form later in a WDM cosmology, they have lower central densities than in CDM. Mass functions are converted from CDM
# to WDM and the completeness corrections calculated perviously are valid for WDM as well.

# Concerning bottom right of figure 7:

# WDM thermal relics with masses from 2 to 12 kev are modeled. There is a big mismatch at low sigma_los. There may be a need to
# include baryonic cores in massive WDM satellites to correct mismatch.

# Concerning figure 8:


# In[ ]:


# 4.4 Results: Summary

# In order to match the observed VF, the following is required:

# 1. The ultrafaint galaxies (𝜎_los < 10 km/s) have dense, cuspy halos, in line with other work on their central densities.
# 2. Tidal stripping must not have reduced the mass of the average ultrafaint galaxy’s dark matter halo by much more than ∼90%.
# 3. The classical dwarf galaxies’ VF (𝜎_los > 10 km/s) is better matched by halos with baryon-driven cores than to cuspy halos,
# although we still overpredict the VF unless subhalos are stripped of ∼ 90% of their dark matter.
# 4. Dwarf satellites must largely survive close to the center of the MW, or else such large completeness corrections apply that
# we have a too many satellites problem.


# In[ ]:


# 5.1 Discussion

# Mismatch is caused by a discrepency at ~10km/s.
# The SMHM relation is required to determine stellar mass of each halo, which is used to calculate sigma_los
# After adoption of SMHM relations, they recover a result that closely matches the more conservative corrected VF values. 
# It is notable that none of the SMHM relations help resolve the too many satellites problem should the corrections 
# accounting for disk stripping be correct.


# In[ ]:


# 5.2 Matching VFs under large completeness corrections

# Larger completeness corrections (associated with greater tidal stripping) predict as many as several 100s to thousands of
# satellites, particularly those with low sigma_los. Variations in Milky Way mass and later reionization could cause this.

# Fiducial reionization redshift z_re = 9.3 could be too low, as higher reionization redshift would allow more halos to grow


# In[ ]:


# 6 Conclusions

# This paper takes observationally based completeness-corrected MW satellite velocity functions and compares them to
# theoretical predicitons.

# The key result is that there is good consistency between fiducial theoretical CDM+baryon velocity function and the most
# conservative completeness-corrected velocity function. 


# In[ ]:


# What to do now:

# reprodude NFW and cored in figure 4

# try to simulate a poisson process

# use mcmc to do compare

# try to answer: 

# 1. How do we evolve from a very homogeneous Universe to our present inhomogeneous state? [A qualitative understanding is enough]
# 2. How do we quantify the amount of structure with size about R?
# 3. How do Dark Matter properties affect the answer to the questions above?


# In[1]:


# 1. After the Big Bang, baryons (electrons and protons at this time) and photons were moving at relativistic speeds inside a  
# small, homogeneous universe that was rapidly expanding. Even though the universe was homogeneous, there were extremely small
# inhomogeneties. Dark matter, which composes 80% of all mass, began to start clustering at this time around the inhomogeneties.
# Eventually, the universe expanded enough that the baryons were able to stop moving at relativistic speeds, and they too began
# to cluster around the already existing gravitational potential wells that were formed by dark matter. Baryons were also formed
# into hydrogen and helium atoms at this point. At this time, the universe was opaque because the non-ionized hydrogen effectively
# absorbed all light and scattered it. As stars formed, they emitted high energy UV radiation that ionized the hydrogen, separating
# protons and electrons once again. This allowed light to shine throughout the universe, but also affected galaxy formation by
# its effect on the cooling and collapse of gas clouds. After the reionization, many stars were able to form, and eventually
# galaxies and galaxy clusters followed. All of these collapses were happening inside dark matter potential wells. Galaxies began to
# collide, forming large elliptical and spiral galaxies. Many dwarf galaxies also began to fall in around these larger galaxies. 
# 2. 
#
#

