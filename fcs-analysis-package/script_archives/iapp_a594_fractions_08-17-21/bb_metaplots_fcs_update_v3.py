#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 21:10:19 2021

@author: Bryan Bogin (bab226)

Purpose: To create nice bar plots from FCS analysis "Batch_fitting_FCS_edit"
This script as well as the Batch fitting scripts are adapted from Greg Gomes.

ARCHIVED: Note this script is archived because it has been designed for a 
specific usage and analysis of FCS data. For a general version, see the main 
directory of the fcs-analysis-package.

"""

#Required Packages
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os 
import pandas as pd
import glob
from uncertainties import ufloat
from uncertainties.umath import *

 
###### Functions ######

def ReadData(filename):
    
    '''
    WARNING: I ADDED A LIMITER (ONLY RETURN FIRST FIVE) BECAUSE I WANTED EACH TO HAVE ONLY 5 DATA POINTS!
    '''
    sample_size = 3   #FIXME
    t = []
    D = []
    err = []
    
    hl=1
    with open(filename, newline = '') as f:                                                                                          
        reader = csv.reader(f, delimiter='\t')
        for j in range(hl):
            #just skip headerlines, already stored
            next(reader)
        for row in reader:
            #Read remaining lines (row at a time)
            
            t.append(float(row[0]))
            D.append(float(row[1]))
            err.append(float(row[2]))


    return np.array(t)[:sample_size], np.array(D)[:sample_size], np.array(err)[:sample_size]

def get_error(n, error):
    '''Error Propagation'''
    error = (np.sqrt(error[0]**2+error[1]**2+error[2]**2))/3
    return float(error)

def average_data(path, files):
    '''Imports data and calculates average using get_error above.
    Saves all data together in a new list.'''
    sample_size = 3   #FIXME
    new_list = []
    
    for i in range(0, len(files)):
        t, D, err = ReadData(path + files[i])
        D = np.array(D)
        Davg = D.mean()
        avg_err = get_error(sample_size, err)
        new_list.append([Davg, avg_err])
    
    return np.array(new_list)

def get_cmap(n, name='viridis'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def Plot_bars(data, label, ylabel, plottitle, savetitle, mode):
    '''Greg's function to that plots bar graphs after loading data.
    Modified with a wrapper, new color scheme, and propagated error bars.'''
    #Plotting Settings
    n = len(data[0])
    color = get_cmap(n+1000)
    
    width = 0.15
    figwidth = 3.42*2
    fig = plt.figure(figsize=(figwidth,figwidth/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    
    rep = list(['Replicate A', 'Replicate B'])
    r1 = np.arange(len(rep))
    
    #Get uncertainty
    
    
    #Plot mean +/- Standard Error of Mean (SEM)
    #d = figwidth/2
    #x = d/r1
    
    for bars in range(0, len(label)):
        disp = 0.15*bars-0.3
        sel1 = bars*2
        sel2 = bars*2+2
        rects1 = ax.bar(r1+disp, data[sel1:sel2,0], color = color(bars*100), yerr = data[sel1:sel2,1]/np.sqrt(n),
               width = width, edgecolor = 'black', label=label[bars])
        ax.bar_label(rects1, padding=1, fmt = '%.2f')
    # rects2 = ax.bar(r1 - 0.1, var2[:,0], color = color(350), yerr = var2[:,1]/np.sqrt(n),
    #               width = width, edgecolor = 'black', label=label[1])
    # rects3 = ax.bar(r1 + 0.1, var3[:,0], color = color(550), yerr = var3[:,1]/np.sqrt(n),
    #         width = width, edgecolor = 'black', label=label[2])
    # rects4 = ax.bar(r1 + 0.3, var4[:,0], color = color(800), yerr = var4[:,1]/np.sqrt(n),
    #         width = width, edgecolor = 'black', label=label[3])
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(ylabel)
    ax.set_title(plottitle)
    ax.set_xticks(r1)
    ax.set_xticklabels(rep)
    
    
   
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    
    #Lit. values for D
    if mode == 0:
        ax.axhline(y=319,color='red',linestyle='--', linewidth=2, label="monomer (lit value)")
        ax.axhline(y=323,color='red',linestyle=':', linewidth=1)
        ax.axhline(y=314,color='red',linestyle=':', linewidth=1)
        
        ax.axhline(y=139,color='black',linestyle='--', linewidth=2, label="oligomer (lit. value)")
        ax.axhline(y=135,color='black',linestyle=':', linewidth=2)
        ax.axhline(y=144,color='black',linestyle=':', linewidth=2)
    
    if mode == 1:
        ax.axhline(y=0.7,color='red',linestyle='--', linewidth=2, label="monomer (lit value)")
        ax.axhline(y=0.65,color='red',linestyle=':', linewidth=1)
        ax.axhline(y=0.75,color='red',linestyle=':', linewidth=1)
        
        ax.axhline(y=1.6,color='black',linestyle='--', linewidth=2, label="oligomer (lit. value)")
        ax.axhline(y=1.55,color='black',linestyle=':', linewidth=2)
        ax.axhline(y=1.65,color='black',linestyle=':', linewidth=2)
    else:
        print("")
   
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    
    
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    #plt.savefig("/Users/bab226/Documents/yale_research/iapp/fcs/fcs-analysis-package/Figures/" + savetitle + ".pdf", dpi=300, transparent=False, bbox_inches='tight')
        
def Get_filenames(path, suffix):
    D_array = []
    print("\n Getting your files... \n")
    for file in glob.glob(path + suffix):   #Ex. "*D1c_DD.dat" to get all dat files with this suffix
        D_array.append(file)
    D_sorted = np.sort(D_array)
    print(D_sorted)


###### Main ######
import os
path = '/Users/bab226/Documents/yale_research/iapp/fcs/fcs-analysis-package/Results/iapp_a594_fractions_08-17-21/'
D1c_files = path+'*D1c_AA.dat'
Rh1c_files = path+'*Rh1c_AA.dat'
C1c_files = path+'*C_1c_AA.dat'

D1c = []
Rh1c = []
C1c = []

#Import files with diffusion constants and Rh values for 1 component model.
for name in glob.glob(D1c_files):
    name = os.path.basename(name)
    D1c = np.append(D1c, name)
    
for name in glob.glob(Rh1c_files):
    name = os.path.basename(name)
    Rh1c = np.append(Rh1c, name)
    
for name in glob.glob(C1c_files):
    name = os.path.basename(name)
    C1c = np.append(C1c, name)

#Sort files
D1c.sort()
Rh1c.sort()
C1c.sort()

print(D1c)
print(Rh1c)
print(C1c)

#Plot D constants
D = average_data(path, D1c)

label = ["a594", "frac08", "frac09", "frac10", "frac11"]
ylabel = 'D (um^2/s)'

Plot_bars(D, label, ylabel, 'Diffusion coefficient', 'diffusion_coefficient', 0)

#Plot Rh
Rh = average_data(path, Rh1c)

label = ["a594", "frac08", "frac09", "frac10", "frac11"]
ylabel = 'Rh (nm)'

Plot_bars(Rh, label, ylabel, 'Rh', 'Rh', 1)

#Plot C
C = average_data(path, C1c)

label = ["a594", "frac08", "frac09", "frac10", "frac11"]
ylabel = 'Conc (nM)'

Plot_bars(C, label, ylabel, 'C', 'C', 3)

##https://www.geeksforgeeks.org/plotting-multiple-bar-charts-using-matplotlib-in-python/?fbclid=IwAR3yOEpVAQz_smHQkzQEfnbuek5AgovowG2db6qcuyhJ1TKVqATCtZkztZE
#how to calculate average error
#move drawing function to separate file



