
# coding: utf-8

# In[51]:


import numpy as np
from scipy import integrate

H0 = 95 #*u.pc #the sacale height of the galaxy
R0 = 2.9e3 #*u.pc # the scale radius
z = 24 #*u.pc #galaxy mid-plane height of the sun
l= 8700 #*u.pc #the placement of the sun 
R_sun = 8700 #*u.pc #we used 8700 bc what was in adams paper not what 8500 what says in the blog
R_super = 3e-2 #/u.yr #galaxy rate of the SN 
R_sgrb=8e-6 #/u.yr
R_lgrb=(1/(1100e6)) #/u.yr
J_kill = 1e5 #(u.J/u.m**2)

class rates:
    
    H0 = 95 #*u.pc #the sacale height of the galaxy
    R0 = 2.9e3 #*u.pc # the scale radius
    z = 24 #*u.pc #galaxy mid-plane height of the sun
    l= 8700 #*u.pc #the placement of the sun 
    R_sun = 8700 #*u.pc #we used 8700 bc what was in adams paper not what 8500 what says in the blog
    R_super = 3e-2 #/u.yr #galaxy rate of the SN 
    R_sgrb=8e-6 #/u.yr
    R_lgrb=(1/(1100e6)) #/u.yr
    J_kill = 1e5 #(u.J/u.m**2)

    
    def Adams_rate(r, R_sn): 
        q0 = (R_sn)/(4*np.pi*H0*R0**2) #normalization 
        q_adams = q0*np.exp(-R_sun/R0)*np.exp(-(np.abs(z))/H0)
        return q_adams*((4/3)*np.pi*(r**3))


    def integral(r, R_sn):
    
        def q_adams(R,z):
            return ((R_sn)/(4*np.pi*H0*R0**2))*np.exp(-R/R0)*np.exp(-(np.abs(z))/H0)
    
    #f = lambda b, l, r: q_adams(R,z)*r**2*np.cos(b) 
    
        def integrand(b,l,r): 
            R = (R_sun**2+(r*np.cos(b))**2-2*R_sun*np.cos(b)*np.cos(l*r))**(1/2)
            z = r*np.sin(b)
            q=q_adams(R,z)
            return q*r**2*np.cos(b) 

    #f_TEST = lambda b,l,r: (r**2)*np.cos(b)
        i = integrate.tplquad(integrand, 0,r, lambda r: 0, lambda r: 2*np.pi,
                                    lambda r, l: -np.pi/2, lambda r, l: np.pi/2)      
        return (i)
    
    def E(r_kill, J_kill):
        return 4*np.pi*r_kill**2*J_kill
    def fluence(E, r):
        return E/(4*np.pi*r**2) #units of J/m^2




# In[87]:

rates.integral(10, R_super)


# In[88]:
'''

integral(20, R_super)


# In[89]:


integral(100, R_super)


# In[90]:


integral(200, R_super)


# In[91]:


integral(400, R_super)


# In[92]:


integral(800, R_super)


# In[93]:


integral(50, R_super)


# In[100]:


integral(3000, R_super)


# In[101]:


integral(1000, R_super)


# In[102]:


8.635488337547774e-05*(9)

'''