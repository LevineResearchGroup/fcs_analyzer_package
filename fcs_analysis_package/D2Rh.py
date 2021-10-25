#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 16:43:54 2021

@author: bab226
"""
import numpy as np

from uncertainties import ufloat
from uncertainties.umath import *


D = ufloat(60,1)
D = D*(1E-6)**2 #now m^2/s 

T = 21.5 + 273.15
A=2.414E-5 # Pa*S
B=247.8 #
C=140

eta=A*10**(B/(T-C))

k = 1.3807E-23

Rh = k*T / (6*np.pi*eta*D)

print(Rh*1E9)

Rh = 1.6E-9

D = k*T/(6*np.pi*eta*Rh)
print(D*1E12)