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

#distribution of radii
r = np.logspace(1,10, num=100, base = 10)
#rate of supernovae
R_sn = rates.R_super
R_lgrb = rates.R_lgrb
R_sgrb = rates.R_sgrb
R = [R_sn, R_lgrb, R_sgrb]


rate_sn = []
rate_lgrb = []
rate_sgrb = []
rate = [rate_sn, rate_lgrb, rate_sgrb]

error_sn = []
error_sgrb = []
error_lgrb = []
error = [error_sn, error_lgrb, error_sgrb]

for r_x in r:
    for R_x in R:
        for i in range(0,2):
            rate[i].append(rates.integral(r_x, R_x)[0])
            error[i].append(rates.integral(r_x, R_x)[1])
        
data = Table([rate_sn, error_sn, rate_lgrb, error_lgrb, rate_sgrb, error_lgrb], names = ['0'])
