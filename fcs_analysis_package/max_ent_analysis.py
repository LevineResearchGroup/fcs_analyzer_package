#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 21:52:01 2021

Purpose: Maximum entropy method analysis. 

For use with MEMFCS software (Quickfit3)
JW. Krieger, J. Langowski (2015): QuickFit 3.0 

The Maximum Entropy approach 
@author: bab226 (Bryan Bogin)

Steps to prepare file:
1. Save autocorrelation curve as a .dat file (use ConvertSPT_to_SimpleTxt.py)
2. Analyze with MaxEnt model.
3. Export data as csv.

What this script does:
4. Imports csv file from Quickfit3.
5. Fits data with gaussian (log-normal distribution)
6. Determines, mean, sd, geometric mean, and geometric SD.
7. Plots fit and data, enables option to compare with conventional fit.

"""

#Load packages
import csv
import numpy as np
import pandas as pd
from numpy import loadtxt
import matplotlib.pyplot as plt
import lmfit
import glob
import os
import re

import FCS_fitfunc as ff

def gauss_fit(log_D, prob, log_D_min, log_D_max):
    #Fit Gaussian using lmfit
    a_guess = np.max(prob)
    mu_guess = np.mean(log_D)
    sigma_guess = np.std(log_D)
    
    print('Guessing a = %f ms' %(a_guess))
    print('Guessing mu = %f ms' %(mu_guess))
    print('Guessing sigma = %f ms' %(sigma_guess))
        
    model = lmfit.Model(ff.gaussian)
    params = model.make_params(a=a_guess, mu=mu_guess, sigma=sigma_guess)
    params['a'].set(min=1E-10, value = a_guess, vary = True)
    params['mu'].set(min=1e-10, value = mu_guess, vary = True)
    params['sigma'].set(min=1e-10, value = sigma_guess, vary = True)
    
    fitres = model.fit(prob[np.logical_and(log_D < log_D_max, log_D > log_D_min)], log_D=log_D[np.logical_and(log_D < log_D_max, log_D > log_D_min)], params=params, method='least_squares',
                           weights=None)
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])

    print(fitres.fit_report())
    return fitres

D_maxent_log = []
D_maxent_log_sigma = []

# measurement_time = float(re.findall(r'T\d+',m['name'])[0][1:])
#Directory and file path
path = "/Users/bab226/Documents/yale_research/iapp/fcs/fcs_analyzer_package/fcs_analysis_package/Data/maxent_quickfit3/"
# name = "test.csv"
for name in glob.glob(path + "40k*"):
    name = os.path.basename(name)[:-4]
    time = float(re.findall(r't\d+', name)[0][1:])
    
    key = "DD"
    #Choose range of diffusion coefficients to select (um^2/s)
    D_min = 10
    D_max = 1000
    
    #Import csv with pandas
    data = pd.read_csv(path+name+".csv")   #read the csv file (put 'r' before the path string to address any special characters in the path, such as '\'). Don't forget to put the file name at the end of the path + ".csv"
    df = pd.DataFrame(data, columns=[data.columns[11], data.columns[12]])    #Select Diffusion coeff. and probability data
    df = df.replace(r'^\s*$', np.NaN, regex=True)   #Needed to replace blanks to nan (will not be able to convert to float if you do not do this)
    
    #Save data to variables
    D = df[df.columns[0]].to_numpy(dtype='float', na_value=np.nan)
    prob = df[df.columns[1]].to_numpy(dtype='float', na_value=np.nan)
    
    #Select specific range of Diffusion coefficients
    mask = np.logical_and(D > D_min, D < D_max)
    
    #Mask to select specific D with high probability
    D = D[mask]
    prob = prob[mask]
    # if len(D) == len(prob):
    #     width = 3.42
    #     fig = plt.figure(figsize=(width,width/1.62))
    #     ax = fig.add_axes([0, 0, 1, 1])
    #     ax.plot(D, prob, linewidth =1, linestyle = '-', marker = 'o', markersize = 1, color = 'k')
    #     # ax.set_xscale('log')
    #     ax.set_xlabel("D (um^2/s)")
    #     ax.set_ylabel("probability")
    #     plt.title('Maximum Entropy Distribution for D')
    
    # else:
    #     print("Huston we have ap problem...")
        
    #Get statistics
    D_mean = np.mean(np.log(D))
    D_median = np.median(np.log(D)) 
    D_std = np.std(np.log(D))
    # D_mean = 1
    # D_median = 2
    # D_std = 0.5
    
    print("""
    mean = %s
    median = %s
    std = %s
          
    """ %(D_mean, D_median, D_std))
    
    #Fit data with gaussian
    # fitres = gauss_fit(np.log(D), prob, np.log(D_min), np.log(D_max))
    fitres = gauss_fit(np.log(D), prob, np.log(D_min), np.log(D_max))
    
    mean = fitres.values['mu']
    std = fitres.values['sigma']
    redchi = fitres.redchi
    
    #Append and save data
    D_maxent_log = np.append(D_maxent_log, mean)
    D_maxent_log_sigma = np.append(D_maxent_log_sigma, std)
    
    #Plot Data and fit
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.plot(D, prob, linewidth =2, linestyle = '-', marker = 'o', markersize = 3, color = 'k')
    ax.set_xscale('log')
    ax.set_xlabel("D (um^2/s)", fontsize=14)
    ax.tick_params(which="major", direction="out",length=6, width=3)
    ax.tick_params(which="minor", direction="out",length=4, width=1.5)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.set_ylabel("Probability", fontsize=14)
    plt.title('MaxEnt Analysis at time %s min' %(time), fontsize=14)
    ax.plot(D, fitres.best_fit, '--', label='best fit')
    
    #Compare to regular D average
    x = 52
    s = 5
    plt.axvspan(x-s, x+s, color='blue', alpha=0.5, edgecolor="black", linewidth=2, label="D1c range")
    ax.legend()
    plt.savefig('./Figures/' + name + '_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    print("""

    ###############
    
    Distribution of D (D +/- geometric SD)
    D = %.3f +/- %.3f
    chi-squared = %.8f
        
    ###############
    
    """ %(np.e**mean, np.e**std, redchi))

D_geom_mean = np.e**D_maxent_log
D_geom_sigma = np.e**D_maxent_log_sigma
    
print("""
    #########  Summary #########
    The geometric mean Diff. coefficient is %.4f +/- %.4f
    The mean geometric sigma is %.4f +/- %.4f.
""" %(np.mean(D_geom_mean), np.std(D_geom_mean), np.mean(D_geom_sigma), np.std(D_geom_sigma)))













