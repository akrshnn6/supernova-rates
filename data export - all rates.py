# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 16:17:32 2019

@author: ashvi
"""

from astropy.io import ascii
import numpy as np
from astropy.table import Table, Column
from astropy import units as u
from integral_3D import rates
#%%
#distribution of radii in pc from 1 (essentially the center) to 2.62^10
r = np.logspace(1, 2.62, num=100, base = 10)
R_sn = rates.R_super
R_lgrb = rates.R_lgrb
R_sgrb = rates.R_sgrb
R = [R_sn, R_lgrb, R_sgrb]

#rates for each explosion
rate_sn = []
rate_lgrb = []
rate_sgrb = []
rate = [rate_sn, rate_lgrb, rate_sgrb]

#integration error for each explosion
error_sn = []
error_sgrb = []
error_lgrb = []
error = [error_sn, error_lgrb, error_sgrb]
#%%

for r_x in r:
    for i in range(0,3):
        rate[i].append(rates.integral(r_x, R[i])[0])
        error[i].append(rates.integral(r_x, R[i])[1])
        
#%%
        

#%%
data_points = []
for i in range(0,3):
    data_points.append(rate[i])
    data_points.append(error[i])
    
name = ['Rate of Supernovae','Error of Supernovae','Rate of LGRB','Error of LGRB','Rate of SGRB','Error of SGRB']
data = Table(data_points, names = name)
data.add_column(Column(r, name = 'Distance'), index = 0)

#%%
ascii.write(data, 'rates.csv', format = 'csv', fast_writer = False, overwrite = True)


