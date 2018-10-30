#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt


# In[2]:


'''
Exact propagation in 2D. 
Pure python loops. First attempt at getting the logic correctly, optimized using cython/numba later.
'''
def exact_prop_2D(in_wave,out_wave,L_in,L_out,wavel,z):
    pi = np.pi
    
    '''
    Build the input and output domains from the input
    '''
    N_in_x = np.shape(in_wave)[0]
    N_in_y = np.shape(in_wave)[1]
    in_domain_x = np.linspace(-L_in/2,L_in/2,N_in_x)
    in_domain_y = np.linspace(-L_in/2,L_in/2,N_in_y)
    
    
    N_out_x = np.shape(out_wave)[0]
    N_out_y = np.shape(out_wave)[1]
    out_domain_x = np.linspace(-L_out/2,L_out/2,N_out_x)
    out_domain_y = np.linspace(-L_out/2,L_out/2,N_out_y)
    
    step_in_x = L_in/N_in_x
    step_in_y = L_in/N_in_y
    
    '''
    Outer loops over i,j -> loop over output array
    Inner loops over p,q -> loop over input array
    For each ouput point, calculate the contribution from each input point and sum 
    '''
    for i in range(N_out_x):
        for j in range(N_out_y):
            x1 = out_domain_x[i]
            y1 = out_domain_y[j]
            for p in range(N_in_x):
                for q in range(N_in_y):
                    x  =  in_domain_x[p]
                    y  =  in_domain_y[q]        
                    f  = in_wave[p][q]
                    out_wave[i][j] += f*np.exp(((-1j*pi)/(wavel*z))*((x-x1)**2+(y-y1)**2))
    '''
    Finally scale the output
    '''
    out_wave *= ((1/np.sqrt(1j*wavel*z))*step_in_x)
    return


# In[3]:


wavel = 0.5*10**(-6)
pi = np.pi
z = 100000
N = 100
L_in  = 5e-1

in_wave = np.zeros((N,N))
in_wave[int(N/2)-int(N/8):int(N/2)+int(N/8),int(N/2)-int(N/8):int(N/2)+int(N/8)] = 1
out_wave_exact = np.zeros((N,N),dtype='complex128')


# In[4]:


exact_prop_2D(in_wave,out_wave_exact,L_in,L_in,wavel,z)


# In[ ]:


f, (ax1,ax2) = plt.subplots(1,2)
ax1.imshow(in_wave,cmap='jet')
ax2.imshow(np.abs(out_wave_exact),cmap='jet')
plt.show()


# In[ ]:


get_ipython().run_line_magic('timeit', '-r 3 exact_prop_2D(in_wave,out_wave_exact,L_in,L_in,wavel,z)')


# In[ ]:




