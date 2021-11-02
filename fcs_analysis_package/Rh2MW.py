# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 14:06:18 2020

@author: gwg24
"""

import numpy as np
import matplotlib.pyplot as plt


MW = np.array([14044, 25665, 42910, 69322, 157368, 239656, 489324, 606444])
R = np.array([1.64, 2.09, 3.05, 3.55, 4.81, 5.2, 6.1, 8.5])

dex_mw = np.array([9.5, 19.5, 39.1, 73, 110, 250, 500, 2000])
dex_rh = np.array([1.86, 3.24, 4.78, 6.49, 7.82, 11.46, 15.90, 26.89])
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1304934/ 
#Dextran size and Rh relation.

m, b = np.polyfit(np.log(dex_mw),np.log(dex_rh),1)

plt.figure()
plt.plot(np.log(dex_mw), np.log(dex_rh), 'o')
plt.plot(np.log(dex_mw), m*np.log(dex_mw) + b , '-')
plt.xlabel('log MW')
plt.ylabel('log Rh')

# plt.yscale('log')
# plt.xscale('log')

MW_0 = 90000 #Da molecular weight
print(np.exp( m*np.log(MW_0) + b))

R_0 = 4.09 #nm stokes radius
MW_est = 1000*np.exp((np.log(R_0)-b)/m)

print(MW_est)
print(MW_est/3900)

# Rmin = 0.066*(MW_0**0.333)
# print(Rmin) #for a perfectly smooth, uniformly packed, spherical -- i.e., not realistic.

