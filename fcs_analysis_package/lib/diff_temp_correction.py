#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 17:04:00 2021

Temp correction for D

@author: bab226
"""

import numpy as np

from uncertainties import ufloat
from uncertainties.umath import *


D_ref = ufloat(374, 40)      

temperature_ref = ufloat(27.85, 0.5)   # temperature at which reference D was taken (celsius)
temperature_lab = ufloat(25, 0.5)    # our labs temeprature (celsius)


T_lab = temperature_lab + 273.15
T_ref = temperature_ref + 273.15
A=2.414E-5 # Pa*S
B=247.8 #
C=140;

eta_lab=A*10**(B/(T_lab-C)) # Pa*s
eta_ref=A*10**(B/(T_ref-C)) # Pa*s

# k = 1.3807E-23 #Boltzman constant J/K
D_ref_lab = D_ref * (T_lab/ eta_lab) * (eta_ref / T_ref) #diffusion coeffecient of reference compoound at actual lab temperature
print(D_ref_lab)
