#!/usr/bin/env python
# coding: utf-8

# In[12]:


# N times throwing a dice experiment

import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

def dicethrow(N):
    data = np.array([])
    for i in range(N):
        samples = np.random.randint(1, 7, 1000)
        count = np.count_nonzero(samples == 3)
        data = np.append(data, count)
    plt.hist(data, bins=16, align='left', rwidth=0.8, color = 'red', edgecolor='black', density='True')
    plt.title('Histogram of Dice Throws')
    plt.xlabel('# of times counted')
    plt.ylabel('# of 3 counts')
dicethrow(100000)


# In[13]:


def poisson(N):
    data = np.random.poisson(1000/6, N)
    plt.hist(data, bins=16, align='left', rwidth=0.8, color = 'red', edgecolor='black', density='True')
    plt.xlabel('# of times counted')
    plt.ylabel('# of 3 counts')
poisson(100000)


# In[ ]:




