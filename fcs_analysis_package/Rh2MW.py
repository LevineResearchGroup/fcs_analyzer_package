# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 14:06:18 2020

@author: gwg24
"""

import numpy as np
import matplotlib.pyplot as plt


MW = np.array([14044, 25665, 42910, 69322, 157368, 239656, 489324, 606444])
R = np.array([1.64, 2.09, 3.05, 3.55, 4.81, 5.2, 6.1, 8.5])

m, b = np.polyfit(np.log(MW),np.log(R),1)

plt.figure()
plt.plot(np.log(MW), np.log(R), 'o')
plt.plot(np.log(MW), m*np.log(MW) + b , '-')
plt.xlabel('log MW')
plt.ylabel('log Rs')

# plt.yscale('log')
# plt.xscale('log')

MW_0 = 90000 #Da molecular weight
print(np.exp( m*np.log(MW_0) + b))

R_0 = 0.72 #nm stokes radius
MW_est = np.exp((np.log(R_0)-b)/m)

print(MW_est)
print(MW_est/3900)

# Rmin = 0.066*(MW_0**0.333)
# print(Rmin) #for a perfectly smooth, uniformly packed, spherical -- i.e., not realistic.

