# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:24:12 2018

@author: Ashvini
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u 

#density vs R for the midplane
#Values found from the paper
#A*np.exp(R/R_d)*np.exp(-np.abs(z))/H
A = 30 #A_v
R_d = 2.9 #kpc
H = ((95*u.pc).to(u.kpc)).value #kpc

#Only accounts for large scales
def density_ccSN(R, z):
    return (A*np.exp(-np.abs(R)/R_d)*np.exp(-np.abs(z))/H)
    

#z = 0 line plot
H_sun =(( 24*u.pc).to(u.kpc)).value #kpc
R_sun = ((8.7)*(u.kpc)).value #kpc
R = np.linspace(-10, 10, num = 10**3)
z = np.linspace(-H, H, num = 10**3)
z_sun = H_sun*u.kpc

plt.figure()
plt.title("Density of ccSN vs Galactic radius (R), z = 0")
plt.ylabel("Density of ccSN")
plt.xlabel("Radius (kpc)")

plt.plot(R, density_ccSN(R, 0),  label ='z = 0')
plt.plot(R, density_ccSN(R, H_sun), label = 'z = Sun Height')
plt.plot(R, density_ccSN(R, H), label = 'z = H')
plt.plot(R, density_ccSN(R, 2*H), label = 'z = 2H')
plt.plot(R_sun, density_ccSN(R_sun, H_sun), 'ro', label = 'Sun/Earth Distance')
plt.legend()

plt.figure()
plt.title("Density of ccSN vs Galactic height (z)")
plt.ylabel("Density of ccSN")
plt.xlabel("Height (kpc)")
plt.plot (z, density_ccSN(0, z), label = 'R = 0')
plt.plot (z, density_ccSN(R_sun, z), label = 'R = R_sun')


plt.title('R-z Plane Contour')
plt.xlabel('Height')
plt.ylabel('Densities')

rho_relation = [1.0, 0.9, 0.75, 0.5, 0.25]
max_density = density_ccSN(R, H)
for value in rho_relation:
    plt.plot(z, value*max_density, label = 'Rho/Rho max= ' +str(value))
    
plt.legend()
plt.legend(bbox_to_anchor = (1, 0.8), loc = 0)

plt.figure()
y = []
for value in z:
    y.append(density_ccSN(R, value))
plt.contour(y)
#plt.ylim(-10, 10)
#plt.xlim(-10, 10)
