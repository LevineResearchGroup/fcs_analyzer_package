# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:41:51 2021

@author: gwg24
"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import math
import os
import subprocess
import scipy as sp
import time



import sys
sys.path.append("./lib") #point python to location containing the below three modules
import FCS_fitfunc as ff
import SPT_reader as spt
import FCS_helpful as fcs


directory_to_bme = './BME-master/'
sys.path.append(directory_to_bme)
import bme_reweight as bme


def diffusion_3d(timelag, tau_diff, A0, Ginf, kappa):
    #1component diffusion 3-D
       return ( Ginf + A0 * 1/(1 + timelag/tau_diff) *
               1/np.sqrt(1 + timelag/(tau_diff*(kappa**2))))

def GaussianDist(x, mu, sig):
    
    return ( (sig*np.sqrt(2*np.pi)**(-1)) * np.exp(-(1/2) * ((x - mu)/2)**2) )



calc_data = 'logtau_prior.dat'

name = 'ADDL_Tube82_aliquot3'
m = spt.Read_FCS('./Data/ADDL_Triplicate_11_05_2021/' + name +'.dat')
# m = spt.Read_FCS(r'C:\Data_LevineLab\FRET-FCS\Data\ADDL_Triplicate_11_05_2021\R110_calibration.dat')
m = m[0]

t_exp = m['DD']['time']
G_exp = m['DD']['G']
err_exp = m['DD']['err']

G_exp = diffusion_3d(t_exp, tau_diff = 0.0055, A0 = np.mean(G_exp[:10]), Ginf = 0, kappa = 5.8)

exp_data = 'BME_format_FCS.dat'
with open(exp_data,'w',newline='\n') as f:
    f.write('# DATA=SAXS PRIOR=GAUSS \n')
    for i in range(len(t_exp)):
        # f.write(labels[i]+'\t')
        f.write("%f \t" % (t_exp[i]))
        f.write("%f \t" % (G_exp[i]))
        f.write("%f \n" % (err_exp[i]))


kappa = 5.8
t_prior = t_exp
tau_prior = np.logspace(-3,3, 100) #ms


G_prior_array = []
for tau in tau_prior:
    
    G_temp = diffusion_3d(t_prior, tau_diff = tau, A0 = np.mean(G_exp[:10]), Ginf = 0, kappa = kappa)
    G_prior_array.append(G_temp)



G_prior_array = np.array(G_prior_array)

print(np.shape(G_prior_array))

mean_G_prior = np.mean(G_prior_array, axis = 0)

plt.figure()
plt.plot(t_exp, mean_G_prior)
plt.plot(t_exp, G_exp)
plt.xscale('log')

print(np.shape(mean_G_prior))

G_calc_prime = np.interp( t_exp, t_prior, mean_G_prior)
print(np.shape(G_calc_prime))


x = G_calc_prime
y = G_exp
XX = np.vstack((x, np.ones_like(x))).T
scaling = np.linalg.lstsq(XX[:, :-1], y, rcond = -1)[0]  # use [0] to just grab the coefs
print(scaling)


#See Bottaro et al., 2018 BME paper. Need to rescale, and have same number of data points for reweighting to work.
BME_calc = G_prior_array*scaling[0]
fun = sp.interpolate.interp1d(t_prior, BME_calc,axis=1, bounds_error=False, fill_value="extrapolate")
BME_calc_prime=fun(t_exp)
print(np.shape(BME_calc_prime))
# np.savetxt(pth+'SAXS_row_framed_norm.dat',BME_calc_prime,delimiter='\t')

with open(calc_data,'w',newline='\n') as f:
    f.write('#Calculated curves \n')
    for i in range(np.shape(BME_calc_prime)[0]):
        f.write(str(tau_prior[i])+'\t')
        np.savetxt(f,[BME_calc_prime[i,:]],delimiter='\t')
        # f.write("\n".join(BME_calc_prime[i,:]))
        
        
        
# thetas = [10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,500,700,1000,5000]
thetas = [0.05, 0.5, 5, 10, 40, 200, 1000, 5000, 10000] 
sele = [10, 20, 50, 100, 5000]

# initialize reweighting class                                                                                                                                
rew = bme.Reweight()

rew.load(exp_data,calc_data)
plt.figure()
for t in thetas:
    # do the minimization
    chi2_before,chi2_after, srel = rew.optimize(theta=t)
    print(chi2_after,np.exp(srel),t)
    w_opt = rew.get_weights()
  
    plt.plot(tau_prior, w_opt, label = str(t))
       
plt.xscale('log')     
plt.legend()

# #Perform the optimization
# chi2_before,chi2_after, srel = rew.optimize(theta = theta)

# # returns the optimized weights

# plt.figure()
# plt.plot(tau_prior, w_opt)
# plt.xscale('log')        