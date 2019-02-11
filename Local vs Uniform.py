# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 12:49:41 2018

@author: ashvi
"""

import numpy as np
from astropy import constants as c
from astropy import units as u 

R_sn = 3*(u.yr**-2)

R_0 = 2.9*u.kpc
H_0 = 95*u.pc

R_g = 15*u.kpc
H_g = 100*u.pc 

q_0 = R_sn/(4*np.pi*H_0*R_0**2)

def gamma_local(q_0, d):
    return (q_0*(4/3)*np.pi*d**3).to(u.yr**-2)
    
print(gamma_local(q_0, 10*u.pc))
    