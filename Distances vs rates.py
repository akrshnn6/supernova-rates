# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 23:40:24 2018

@author: Ashvini
"""
import numpy as np
from astropy import constants as c
from astropy import units as u
#import matplotlib.pyplot as plt

A_i = 60 #atomic number of iron
M_p = c.m_p #mass of a proton
Xi = 10**-5 #fraction of mass
A = 30 #A_v
R_d = 2.9 #kpc

H = (95*u.pc) #kpc
H_sun =(( 24*u.pc).to(u.kpc)) #kpc

R_sun = ((8.7)*(u.kpc)).to(u.pc) #kpc
R = np.linspace(-10, 10, num = 10**3)

z = np.linspace(-H, H, num = 10**3)
z_sun = H_sun*u.kpc

N_dot = 3e-2
m_u = (58.9348755 * u.u).to(u.kg)

#^60F surface density based on our derivation 
def density_60F(M_eji, D):
    return (M_eji/(16* np.pi* D**2 *A_i*m_u))

#Simple disk approximation
def sn_rate(D, N_dot, R, H):
    return (2*D**3*N_dot)/(3*R**2*H)

N_i = (1.41e6)*u.cm**-2
M_eji =  (10*u.M_sun).to(u.kg)
D = 10*u.pc

m_ejected = [3.6*10**-5, 6.6*10**-5, 1.1*10**-4, 3.6*10**-5, 2.5*10**-5, 1.5*10**-4]*u.M_sun
m_ccSN = [9, 15, 19, 20, 21, 25]*u.Msun
colors = ['red', 'blue' 'green', 'orange', 'purple']

#Radioactivity distance based on our derivation 
def rad_distance(M_eji, N_i):

    return (1/4)*np.sqrt(M_eji/(np.pi * N_i * A_i * m_u))


#%%
densities = []
distances = []
for value in range(len(m_ejected)):
    dens = (density_60F(m_ejected[value], D).to(u.pc**-2))
    densities.append(dens.value)
    dist = (rad_distance(m_ejected[value], (1.41e6)*u.cm**-2).to(u.pc))
    distances.append(dist.value)

#%%
dist_paper = [70, 90, 120, 72, 60, 130]
rates = []
prates =[]
for value in range(len(distances)):
    rate = distances[value]
    prate = dist_paper[value]
    print('Distance = '+ str(distances[value]))
    print(' SN rate per year: '+str(sn_rate(rate, N_dot, R_sun, H_sun)))
    print(str(sn_rate(rate, N_dot, 10**4*u.pc, H_sun)*3)+ ' Supernovae for 3Myr')
    print('')
    print('Distance = '+ str(dist_paper[value]))
    print(' SN rate per year: '+str(sn_rate(prate, N_dot, R_sun, H_sun)))
    print(str(sn_rate(prate, N_dot, 10**4*u.pc, H_sun)*3)+ ' Supernovae for 3Myr')
    print('')
    rates.append(rate)
    prates.append(prates)
