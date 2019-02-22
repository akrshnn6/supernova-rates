
# coding: utf-8

# In[62]:


import numpy as np
#from astropy import units as u
import matplotlib.pyplot as plt 
#import pylab
#import scipy.integrate as integrate
#import scipy.special as special
#import sympy as sp
from integral_3D import rates
'''
# In[53]:


print("Ionizing SN Energy from the table is 5.46*10^32")
      
0.32*10e31+0.22*2e33+0.02*10e32+0.44e32

print ("Ionizing SN Energy from the plot is 8.36*10^35")


# In[54]:


(0.32*10e31+0.22*2e33+0.02*10e32+0.44e32)


# In[55]:


H0 = 95#*u.pc #the sacale height of the galaxy
R0 = 2.9e3#*u.pc # the scale radius
z = 24#*u.pc #galaxy mid-plane height of the sun
R=8700*u.pc #the placement of the sun 

E_sn = #5.46e32*u.J
E_sgrb= #10e43*u.J
E_lgrb= #5e45*u.J


R_super = #3e-2/u.yr #galaxy rate of the SN 
R_sgrb=#8e-6/u.yr
R_lgrb=#(1/(1100e6))/u.yr

R_grb=R_super*10**(-5)/u.yr
'''

# In[56]:

class double_e:
    def Adams_rate(r, R_sn): 
        q0 = (R_sn)/(4*np.pi*rates.H0*rates.R0**2) #normalization 
        q_adams = q0*np.exp(-rates.R_sun/rates.R0)*np.exp(-(np.abs(rates.z))/rates.H0)
        return q_adams*((4/3)*np.pi*(r**3))


# In[57]:
'''

def Adams_rate_int (r, R_sn):
    q0 = (R_sn)/(4*np.pi*H0*R0**2)
    q_adams = q0*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)
    f1 = (np.pi*q_adams*H0*np.exp(-r/H0))
    f2 = (r**2)-(np.exp(r/H0)-1)*(2*H0**2-r**2)+2*r*H0
    return f1*f2


# In[58]:


Adams_rate_int (10*u.pc, R_super)


# In[59]:


Adams_rate (10*u.pc, R_super)


# In[26]:


def F (r, E):
    return  (E/(4*np.pi*r**2))

'''
# In[28]:
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
#E_sn = 1.14e40*u.J
#=10^47 erg
#%%
r = np.logspace(1, 5.477, num=100, base = 10)    
rate_SN = []
rate_LGRB = []
rate_SGRB = []
for val in r:
    rate_SN.append(double_e.Adams_rate(val, rates.R_super))

plt.plot(r, rate_SN)

'''#%%    
SN = []
sgrb = []
lgrb = []
#calculate r from fluence array and use below
for value in r: 
    #f_sn.append(F(value, E_sn).value)
    #f_sgrb.append(F(value, E_sgrb).value)
    #f_lgrb.append(F(value, E_lgrb).value)
    SN.append(Adams_rate(value, R_super).value)
    sgrb.append(Adams_rate(value, R_sgrb).value)
    lgrb.append(Adams_rate(value, R_lgrb).value)

#pylab.plot(4.84e-10,6.08E35, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, SN, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, sgrb, color = 'orange', label = 'SGRB')
ax1.plot(f, lgrb, color = 'red', label = 'LGRB')
ax1.axvline (x=9.523809523809e37,color='green', linestyle='dashed', label = 'Kill Fluence')

#plt.axhline(y=4.8403911e-10, color='cornflowerblue', linestyle='-')

#ax1.set_xlim(10e34,10e50)
#plt.xlim(10,10e50)

#ax1.set_ylim(10e-30,10)
ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_xlabel ('Fluence (J/pc^2)')
ax1.set_ylabel ('Rate (yr)')
ax1.set_title("Fluence versus Rate (Adams Model)", y = 1.15, fontsize = 15)

ax1.legend()

new_ticks = []
for x in ax1.get_xticks():
    new_ticks.append(x*1e-33)
x = []
for value in new_ticks:
    x.append('%1.e'% value)
    
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(x)
ax2.set_xlabel("(J/m^2)")
plt.yscale('log')
pylab.show()


# In[29]:


#r= np.linspace(1,10e60,10000000)
r= np.linspace(10e-16,10e60,1000)

#pylab.plot(4.84e-10,6.08E35, color = 'cornflowerblue', label = 'Supernovae')
plt.plot(F(r,8.36e35),Adams_rate(r,R_super), color = 'cornflowerblue', label = 'Supernovae')
plt.plot(F(r,E_sgrb), Adams_rate(r,R_sgrb), color = 'orange', label = 'SGRB')
plt.plot(F(r,E_lgrb), Adams_rate(r,R_lgrb), color = 'red', label = 'LGRB')
plt.axvline (x=9.523809523809e37,color='green', linestyle='dashed', label = 'Kill Fluence')

#plt.axhline(y=4.8403911e-10, color='cornflowerblue', linestyle='-')

plt.xlim(10e34,10e50)
#plt.xlim(10,10e50)

plt.ylim(10e-30,10)

plt.xlabel ('Log Fluence (J/pc^2)')
plt.ylabel ('Log Rate (yr)')
plt.title("Fluence versus Rate (Adams Model)", fontsize = 15)

plt.xscale('log')
plt.yscale('log')
plt.legend()
pylab.show()

'''