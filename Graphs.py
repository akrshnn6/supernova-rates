# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 20:29:00 2019

@author: ashvi
"""

from astropy.io import ascii 
import matplotlib.pyplot as plt 
from astropy import units as u
from integral_3D import rates
#%%
rates_exp = ascii.read('rates.csv')
r = rates_exp['Distance']

plt.figure(figsize = (10,15))

plt.subplot(311)
plt.plot(rates_exp['Distance'], rates_exp['Rate of Supernovae'], color = 'b')
plt.title('Rate of Supernovae vs Distance')
plt.ylabel('Rate (SN/yr)')
plt.xlabel('Distances (pc)')

plt.subplot(312)
plt.plot(rates_exp['Distance'], rates_exp['Rate of LGRB'], color = 'g')
plt.title('Rate of LGRB vs Distance')
plt.ylabel('Rate (LGRB/yr)')
plt.xlabel('Distances (pc)')

plt.subplot(313)
plt.plot(rates_exp['Distance'], rates_exp['Rate of SGRB'], color = 'r')
plt.title('Rate of SGRB vs Distance')
plt.ylabel('Rate (SGRB/yr)')
plt.xlabel('Distances (pc)')

#overplot uniform model & double exponential model 

plt.tight_layout()
plt.savefig('rates vs distances.png')
#%%
fluences = ascii.read('fluences.csv')
r = rates_exp['Distance']

plt.figure(figsize = (10,15))

plt.subplot(311)
plt.loglog(fluences['Fluence of Supernovae'], rates_exp['Rate of Supernovae'], color = 'b')
plt.title('Fluence of Supernovae vs Distance')
plt.xlabel('Fluence (erg/pc^2)')
plt.ylabel('Distances (pc)')

plt.subplot(312)
plt.loglog(fluences['Fluence of LGRB'], rates_exp['Rate of LGRB'], color = 'g')
plt.title('Fluence of LGRB vs Distance')
plt.xlabel('Fluence (erg/pc^2)')
plt.ylabel('Distances (pc)')

plt.subplot(313)
plt.loglog(fluences['Fluence of SGRB'], rates_exp['Rate of SGRB'], color = 'r')
plt.title('Fluence of SGRB vs Distance')
plt.xlabel('Fluence (erg/pc^2)')
plt.ylabel('Distances (pc)')

plt.tight_layout()
plt.savefig('fluences vs rates.png')

