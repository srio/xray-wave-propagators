#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/sajid/packages/xray-wave-propagators/prop/cython/binaries")
import prop2d_cython_2_loop


# In[2]:


wavel = 0.5*10**(-6)
z = 100000.0
N = 100
L_in  = 0.5

in_wave = np.zeros((N,N))
in_wave[int(N/2)-int(N/8):int(N/2)+int(N/8),int(N/2)-int(N/8):int(N/2)+int(N/8)] = 1
in_wave = np.array(in_wave,dtype=np.complex128)

out_wave = np.zeros((N,N),dtype=np.complex128)

sampling = L_in/N
critical = (wavel*z/L_in)
if sampling>critical:
    print('Use TF')
else :
    print('Use IR')
print('Fresnel Number :', (L_in**2)/(wavel*z))


# In[3]:


prop2d_cython_2_loop.exact_prop_2D_cython(in_wave,out_wave,L_in,L_in,0,0,wavel,z,5)


# In[4]:


f, (ax1,ax2) = plt.subplots(1,2)
ax1.imshow(np.abs(in_wave),cmap='jet')
ax2.imshow(np.abs(out_wave),cmap='jet')
plt.show()


# In[5]:


get_ipython().run_line_magic('timeit', 'prop2d_cython_2_loop.exact_prop_2D_cython(in_wave,out_wave,L_in,L_in,0,0,wavel,z,5)')


# In[ ]:




