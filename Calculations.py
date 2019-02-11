# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:04:40 2018

@author: Ashvini
"""

import numpy as np
from astropy import constants as c
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.integrate import dblquad

A_i = 60 #atomic number of iron
M_p = c.m_p
Xi = 10**-5 #fraction of mass
A = 30 #A_v
R_d = 2.9 *u.kpc #kpc


H = ((95*u.pc).to(u.kpc)) #kpc
N_dot = 3e-2 *u.yr**-1

H_sun =(( 24*u.pc).to(u.kpc)).value #kpc
R_sun = ((8.7)*(u.kpc)) #kpc
R = np.linspace(-10, 10, num = 10**3)
z = np.linspace(-H, H, num = 10**3)
z_sun = H_sun*u.kpc

#^60F surface density based on our derivation 
def density_60F(M_eji, D):
    return (M_eji/(16* np.pi* D**2 *A_i*M_p))

#Simple disk approximation
def sn_rate(D, N_dot, R, H):
    return (2*D**3*N_dot)/(3*R**2*H)

#Volume of a sphere
def vol_sphere(radius):
    return (4/3)*(np.pi)*(radius^3) *u.pc**3

#Only accounts for large scales
def density_ccSN(R, z):
    return (A*np.exp(-np.abs(R)/R_d)*np.exp(-np.abs(z))/H)
    
M_eji =  (10*u.M_sun).to(u.kg)
D = 10*u.pc

m_ejected = [3.6*10**-5, 6.6*10**-5, 1.1*10**-4, 3.6*10**-5, 2.5*10**-5, 1.5*10**-4]*u.M_sun
m_ccSN = [9, 15, 19, 20, 21, 25]
colors = ['red', 'blue' 'green', 'orange', 'purple']

#Radioactivity distance based on our derivation 
def rad_distance(M_ej, N_i):
    return np.sqrt(M_ej.to(u.kg)/(N_i* 16* np.pi * A_i* M_p))


#%%
        
def q_SN(rate, volume):
    return rate/volume 
    
def q_SN_cluster(rate, volume, number):
    return q_SN(rate, volume) / number 

def q_0(N_dot, R_0, H_0):
    if type(R_0) == u.quantity.Quantity:
        R_0 = R_0.value
    if type(H_0) == u.quantity.Quantity:
        H_0 = H_0.value
    a = (dblquad(lambda r,z: np.exp(-r/R_0)*np.exp(-np.abs(z)/H_0), 0, np.inf, lambda z: -np.inf, lambda z: np.inf))
    return (N_dot/(np.pi*a[1]*u.pc**3)).to(u.Myr**-1 * u.pc**-3)
    
R_d = 2.9 *u.kpc #kpc
H = ((95*u.pc).to(u.kpc)) #kpc
N_dot = 3e-2 *u.yr**-1

r = 10
rate = sn_rate(D, N_dot, R_d, H) 
volume= vol_sphere(r)
number = 10
#units of pc^-3*Myr^-1
vol_t = u.pc**-3*u.Myr**-1

print('Homogeneous q_SN for D =' + str(r) + ' pc: ' + str(q_SN(rate, volume).to(vol_t)))
#print('Homogeneous Γ(D) for D =' + str(r) + 'pc, and '+ str(number) + ' supernovae in a cluster: ' + str(q_SN_cluster(rate, volume, number).to(vol_t)))

#print('q_0 = '+ str(q_0(N_dot, R_d, H)))
 #%%
#Supernova density for ejected mass of 10 solar masses w/ kill distance @ 10pc
print('Density: '+str(density_60F(M_eji*Xi, D).to(u.cm**-2)))
#Using given value of 10^13 m^-2 for Ni
Ni = 10**13*u.m**-2
print('Distance: '+str(rad_distance((10**-5*u.M_sun), Ni).to(u.pc)))
print('')

densities = []
distances = []
for value in range(len(m_ejected)):
    dens = (density_60F(m_ejected[value], D).to(u.pc**-2))
    print('Supernova Mass '+str(m_ccSN[value]) + 'M_sun')
    print('Density: '+str(dens))
    densities.append(dens.value)
    dist = (rad_distance(m_ejected[value], (1.41e6)*u.cm**-2).to(u.pc))
    print('Distance: '+ str(dist))
    distances.append(dist.value)
    print('')

#%%
plt.figure()    
plt.scatter(m_ccSN, distances, label = 'Distances')
rates = []  
#plt.scatter(m_ccSN, rates, label = 'Rates')
plt.title('Supernova Mass vs Distances')
plt.xlabel('Supernova Mass ($M_{\odot}$)')

plt.ylabel('Distances(pc)')
dist = np.round(distances,3)
smallest = distances[0]
largest = 0
for value in range(len(m_ccSN)):
    d = distances[value]
    if value>0 and m_ccSN[value]!=20:
        plt.text(m_ccSN[value]+1, d+distances[0]*(-0.005), str(m_ccSN[value]) +  '$M_{\odot}$ $CCSN^a$'  )
        plt.text(m_ccSN[value]+1, d+distances[0]*(-0.1), '$d_{pc}$ = '+str(dist[value])+'pc')
    elif m_ccSN[value] ==20:
        plt.text(m_ccSN[value]+1, d+distances[0]*(0.1), str(m_ccSN[value]) +  '$M_{\odot}$ $CCSN^a$'  )
        plt.text(m_ccSN[value]+1, d+distances[0]*(0.005), '$d_{pc}$ = '+str(dist[value])+'pc')
    else:
        plt.text(m_ccSN[value] + 1, d+distances[0]*(-0.005), str(m_ccSN[value]) +  '$M_{\odot}$ $ECSN^b$'  )
        plt.text(m_ccSN[value] + 1, d+distances[0]*(-0.1), '$d_{pc}$ = '+str(dist[value])+'pc')
    if d<smallest:
        smallest = d
    if d>largest:
        largest = d
plt.ylim(smallest-0.2*smallest, largest+0.2*largest)

#%%
plt.figure()    
plt.scatter(m_ccSN, densities)
plt.title('Supernova Mass vs Density')
plt.xlabel('Supernova Mass ($M_{\odot}$)')

plt.ylabel('Densities ($atoms ^{60}F/pc^{-2}$)')
for value in range(len(m_ccSN)):
    if value>0 and m_ccSN[value]!=20:
#        plt.text(m_ccSN[value] - 7, densities[value]-10, str(m_ccSN[value]) +  '$M_{\odot}$ $CCSN^a$'  )
        plt.text(m_ccSN[value] - 7, densities[value]-22**20, '$F_{obs,60}$ = '+str('%.3e' % densities[value])+'($atoms ^{60}F/pc^{-2}$)')
    elif m_ccSN[value] ==20:
#        plt.text(m_ccSN[value], densities[value]+22, str(m_ccSN[value]) +  '$M_{\odot}$ $CCSN^a$'  )
        plt.text(m_ccSN[value], densities[value]+10**20, '$F_{obs,60}$ = '+str('%.3e' % densities[value])+'($atoms ^{60}F/pc^{-2}$)')
    else:
#        plt.text(m_ccSN[value] - 7, densities[value]-10, str(m_ccSN[value]) +  '$M_{\odot}$ $ECSN^b$'  )
        plt.text(m_ccSN[value] - 7, densities[value]-22**20, '$F_{obs,60}$ = '+str('%.3e' % densities[value])+'($atoms ^{60}F/pc^{-2}$)')
