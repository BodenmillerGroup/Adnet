
# coding: utf-8

# # Benchmark bp-R2 vs traditional measures

# In the paper it has been asked how well bp-R2 compares to more traditional correlation markers.

# In[3]:

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')
get_ipython().magic('load_ext autoreload')


# In[ ]:




# In[4]:

get_ipython().magic('autoreload 1')
get_ipython().magic('aimport library')


# In[ ]:




# In[5]:

n=500


# In[6]:

y, x = library.simulate_perfect_distribution(10000,'independent')
plt.scatter(x, y, alpha=0.1)


# In[7]:

x,y = library.simulate_perfect_distribution(100,'linear')
plt.scatter(x,y)


# In[8]:

x,y = library.simulate_perfect_distribution(100,'square')
plt.scatter(x,y)


# In[9]:

y,x = library.simulate_perfect_distribution(100,'exponential')
plt.scatter(x,y)


# In[10]:

y,x = library.simulate_perfect_distribution(100,'quadratic')
plt.scatter(x,y)


# In[11]:

y,x = library.simulate_perfect_distribution(100,'sine')
plt.scatter(x,y)


# In[12]:


x,y = library.simulate_perfect_distribution(10000,'circumference')
plt.scatter(x,y)


# In[13]:

reload(library)


# In[14]:

x,y = library.simulate_perfect_distribution(100,'square')
plt.scatter(x,y)


# In[15]:

x,y = library.simulate_perfect_distribution(100,'cross')
plt.scatter(x,y)


# In[16]:

x,y = library.simulate_perfect_distribution(100)
plt.scatter(x,y)


# In[ ]:



