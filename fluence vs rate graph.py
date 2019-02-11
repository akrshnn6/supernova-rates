# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 13:03:58 2018

@author: ashvi
"""


#%%
import numpy as np
from astropy import units as u
from astropy.table import Table
#from astropy.table.column import Column
import scipy.special as sps
import matplotlib.pyplot as plt
import pylab
#%%

H0 = 95*u.pc 
R0_adams = 2.9e3*u.pc
R_sn = 0.03*u.yr**-1
R_sun = 8700 * u.pc
z = 24*u.pc

#%%

alpha = 1.09
beta =  3.87
R0_green = R_sun/beta
#%%
#comparision of all models within this cell


a = R_sn/(4*np.pi*sps.gamma(alpha+2)*R0_green**2*H0)
b = (R_sun/R0_green)**alpha
c = np.exp(-R_sun/R0_adams)
d = np.exp(-np.abs(z)/H0)


Hg, Rg = 100*u.pc, 15000*u.pc
q_uniform = R_sn/(2*np.pi*Hg*Rg**2)
q_uni = q_uniform
q_adams = (R_sn)/(4*np.pi*H0*R0_adams**2)*c*d
q_ad = q_adams
c = np.exp(-R_sun/R0_green)
q_green = a*b*c*d
q_gr = q_green
q_s = [q_uniform, q_adams, q_green]    

#returns the rate in units of SN/Gyr
#enter d in units of pc 
def rate(q, d):
    return q*(4/3)*np.pi*(d)**3

#adams dist 175.131 pc
#enter rate in units of yr^-1
sn_rate = 2.6e-6*u.yr**-1
def mean_dist(q, rate):
    return ((3*rate)/(4*np.pi*q))**(1/3)
 
#adams: 2.066e9 yrs
#enter d in units of pc
def mean_time(q, d):
    return 1/(q*(4/3)*np.pi*(d)**3)
#%%
#Killer GRBs
#short bursts from neutron stars merging, long from supernovae
d = 2000*u.pc
#factor of 10**-5 comes from rate of grb/rate of sn find that number 
q_GRB = q_adams*10**-5
rate_SGRB = rate(q_GRB, 2000*u.pc)
rate_GRB = rate(q_GRB, 200*u.pc)

#%%

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

J_kill = 1e5*(u.J/u.m**2)
def E(r_kill):
    return 4*np.pi*r_kill**2*J_kill
def fluence(E, r):
    return E/(4*np.pi*r**2)

SN_kill = 10*u.pc
SGRB_kill = 200*u.pc
LGRB_kill = 2000*u.pc

E_sn = E(SN_kill).to(u.J)
E_sgrb = E(SGRB_kill).to(u.J)
E_lgrb = E(LGRB_kill).to(u.J)

def rate_F(q, E, f):
   return q*(4/3)*np.pi*(E/(4*np.pi*f))**(3/2)
    
rate_sn = []
rate_sgrb = []
rate_lgrb = []

f = np.logspace(0, 10, num=1000, base=10)
for value in f:
    rate_sn.append(rate_F(q_adams, E_sn, value*(u.J/u.m**2)).to(1/u.yr).value)
    rate_sgrb.append(rate_F(q_GRB, E_sgrb, value*(u.J/u.m**2)).to(1/u.yr).value)
    rate_lgrb.append(rate_F(q_GRB, E_lgrb, value*(u.J/u.m**2)).to(1/u.yr).value)
    

#pylab.plot(4.84e-10,6.08E35, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, rate_sn, color = 'cornflowerblue', label = 'Supernovae')
ax1.plot(f, rate_sgrb, color = 'orange', label = 'SGRB')
ax1.plot(f, rate_lgrb, color = 'red', label = 'LGRB')
ax1.axvline (x=1e5,color='green', linestyle='dashed', label = 'Kill Fluence')

#plt.axhline(y=4.8403911e-10, color='cornflowerblue', linestyle='-')

#ax1.set_xlim(0,1e43)
#ax1.set_xlim(0,1e10)

#ax1.set_ylim(10e-30,10)
ax1.set_xscale('log')
ax2.set_xscale('log')
ax1.set_xlabel ('Fluence (J/m^2)')
ax1.set_ylabel ('Rate (yr)')
ax1.set_title("Fluence versus Rate (Adams Model)", y = 1.15, fontsize = 15)

ax1.legend()

new_ticks = []
for x in ax1.get_xticks():
    new_ticks.append(x*9.521e32)
x = []
for value in new_ticks:
    x.append('%1.e'% value)
    
ax2.set_xticks(ax1.get_xticks())
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(x)
ax2.set_xlabel("J/pc^2")
plt.yscale('log')
pylab.show()
#%%



