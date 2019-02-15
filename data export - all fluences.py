# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 14:42:07 2019

@author: ashvi
"""
from astropy.io import ascii
import numpy as np
from astropy.table import Table, Column
from astropy import units as u
from integral_3D import rates

file = ascii.read('rates.csv')

#distribution of radii in pc from 1 (essentially the center) to 2.62^10
r = np.logspace(1, 4.177, num=100, base = 10)

R_sn = rates.R_super
R_lgrb = rates.R_lgrb
R_sgrb = rates.R_sgrb
R = [R_sn, R_lgrb, R_sgrb]

J_kill = (1e5*(u.J/u.m**2)).to(u.erg/u.pc**2)

f_sn = []
f_lgrb = []
f_sgrb = []
f = [f_sn, f_lgrb, f_sgrb]
for rad in r:
    for i in range(0,3):
        f[i].append(rates.fluence(rates.E(R[i], J_kill).value,rad))
        
data_points = f
name = ['Fluence of Supernovae', 'Fluence of LGRB', 'Fluence of SGRB']
data = Table(data_points, names = name)
data.add_column(Column(r, name = 'Distance'), index = 0)
ascii.write(data, 'fluences.csv', format = 'csv', fast_writer = False, overwrite = True)
