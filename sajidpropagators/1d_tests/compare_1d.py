#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import numexpr as ne

from sajidpropagators.prop.propagators_1d import exact_prop, propIR
# get_ipython().run_line_magic('run', '../prop/propagators_1d.py')


# In[2]:

energy = 10000
wavel = (1240/energy)*10**(-9)
#wavel = 0.5*10**(-6)
pi = np.pi
z = 1000e-6
N = 1000
L_in  = 10e-6

in_domain_exact  = np.linspace(-L_in/2,L_in/2,N)
in_wave = np.zeros(N)
in_wave[int(N/2)-int(N/8):int(N/2)+int(N/8)] = 1


# In[3]:


sampling = in_domain_exact[1] - in_domain_exact[0]
critical = (wavel*z/L_in)
print(sampling>critical)
print('Fresnel Number :', (L_in**2)/(wavel*z))


# In[4]:


out_,L_ = propIR(in_wave,L_in/N,L_in,wavel,z)
out_domain_ = np.linspace(-L_/2,L_/2,N)


# In[5]:


N = 1000
in_domain_exact  = np.linspace(-L_in/2,L_in/2,N)
in_wave = np.zeros(N)
in_wave[int(N/2)-int(N/8):int(N/2)+int(N/8)] = 1
out_wave_exact = np.zeros((N),dtype='complex128')
exact_prop(in_wave,out_wave_exact,L_in,L_,wavel,z)


# In[6]:


f, (ax1,ax2,ax3) = plt.subplots(1,3)
ax1.plot(in_domain_exact*1e6,np.abs(in_wave))
ax1.set_xlabel('co-ordinates in um',fontsize = 15)
ax1.set_title('Input', fontsize = 15)
ax2.plot(np.abs(out_wave_exact),'b')
ax2.set_xlabel('co-ordinates in um',fontsize = 15)
ax2.set_title('Output Exact', fontsize = 15)
ax3.plot( out_domain_*1e6,np.abs(out_),'g')
ax3.set_xlabel('co-ordinates in um',fontsize = 15)
ax3.set_title('Output IR', fontsize = 15)
f.set_size_inches(20, 10, forward=True)
f.suptitle('Fresnel Number : '+str((L_in**2)/(wavel*z)),fontsize = 25)
plt.show()


# In[7]:


'''
f, (ax1,ax2) = plt.subplots(1,2)
ax1.scatter(np.abs(out_),np.abs(out_wave_exact))
ax1.set_title('Amplitude overlap', fontsize = 15)
ax2.scatter(np.unwrap(np.angle(out_)),np.unwrap(np.angle(out_wave_exact)))
ax1.set_title('Phase overlap', fontsize = 15)
f.set_size_inches(20, 10, forward=True)
plt.show()
'''






