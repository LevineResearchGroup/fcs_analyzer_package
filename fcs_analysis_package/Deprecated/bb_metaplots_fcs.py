#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 21:10:19 2021

@author: Bryan Bogin (bab226)

Purpose: To create nice bar plots from FCS analysis "Batch_fitting_FCS_edit"
This script as well as the Batch fitting scripts are adapted from Greg Gomes.

"""

#Required Packages
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os 
import pandas as pd
import glob
# from uncertainties import ufloat
# from uncertainties.umath import *
 
###### Functions ######

def ReadData(filename):
    
    '''
    WARNING: I ADDED A LIMITER (ONLY RETURN FIRST FIVE) BECAUSE I WANTED EACH TO HAVE ONLY 5 DATA POINTS!
    '''
    sample_size = 3
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

def Import_data(path, file1, file2):
    '''Imports data using the function above'''
    new_array = []
    t, D, err = ReadData(path + file1)
    new_array.append([np.mean(D), np.std(D,ddof = 1)])
    t, D, err = ReadData(path + file2)
    new_array.append([np.mean(D), np.std(D,ddof = 1)])
    new_array = np.array(new_array)
    return new_array

def get_cmap(n, name='viridis'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def Plot_bars(var1, var2, var3, var4, label, ylabel, plottitle, savetitle, mode):
    '''Greg's function to that plots bar graphs after loading data'''
    #Plotting Settings
    n = len(var1[0])
    color = get_cmap(n+1000)
    
    width = 0.2
    figwidth = 3.42*2
    fig = plt.figure(figsize=(figwidth,figwidth/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    
    labels = list(['Replicate A', 'Replicate B'])
    r1 = np.arange(len(labels))
    
    
    #Plot mean +/- Standard Error of Mean (SEM)
    rects1 = ax.bar(r1 - 0.3, var1[:,0], color = color(50), yerr = var1[:,1]/np.sqrt(n),
              width = width, edgecolor = 'black', label=label[0])
    rects2 = ax.bar(r1 - 0.1, var2[:,0], color = color(350), yerr = var2[:,1]/np.sqrt(n),
              width = width, edgecolor = 'black', label=label[1])
    rects3 = ax.bar(r1 + 0.1, var3[:,0], color = color(550), yerr = var3[:,1]/np.sqrt(n),
              width = width, edgecolor = 'black', label=label[2])
    rects4 = ax.bar(r1 + 0.3, var4[:,0], color = color(800), yerr = var4[:,1]/np.sqrt(n),
              width = width, edgecolor = 'black', label=label[3])
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(ylabel)
    ax.set_title(plottitle)
    ax.set_xticks(r1)
    ax.set_xticklabels(labels)
    
    ax.bar_label(rects1, padding=1, fmt = '%.2f')
    ax.bar_label(rects2, padding=1, fmt = '%.2f')
    ax.bar_label(rects3, padding=1, fmt = '%.2f')
    ax.bar_label(rects4, padding=1, fmt = '%.2f')
   
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    
    #Lit. values for D
    if mode == 0:
        ax.axhline(y=319,color='grey',linestyle='--', linewidth=2, label="monomer (lit value)")
        ax.axhline(y=323,color='grey',linestyle=':', linewidth=1)
        ax.axhline(y=314,color='grey',linestyle=':', linewidth=1)
        
        ax.axhline(y=139,color='black',linestyle='--', linewidth=2, label="oligomer (lit. value)")
        ax.axhline(y=135,color='black',linestyle=':', linewidth=2)
        ax.axhline(y=144,color='black',linestyle=':', linewidth=2)
    
    if mode == 1:
        ax.axhline(y=0.7,color='grey',linestyle='--', linewidth=2, label="monomer (lit value)")
        ax.axhline(y=0.65,color='grey',linestyle=':', linewidth=1)
        ax.axhline(y=0.75,color='grey',linestyle=':', linewidth=1)
        
        ax.axhline(y=1.6,color='black',linestyle='--', linewidth=2, label="oligomer (lit. value)")
        ax.axhline(y=1.55,color='black',linestyle=':', linewidth=2)
        ax.axhline(y=1.65,color='black',linestyle=':', linewidth=2)
    
    ax.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    
    
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    plt.savefig("./Figures/" + savetitle + ".pdf", dpi=300, transparent=False, bbox_inches='tight')
        
def Get_filenames(path, suffix):
    D_array = []
    print("\n Getting your files... \n")
    for file in glob.glob(path + suffix):   #Ex. "*D1c_DD.dat" to get all dat files with this suffix
        D_array.append(file)
    D_sorted = np.sort(D_array)
    print(D_sorted)


###### Main ######
import os
path = './Results/'
D1c_files = './Results/*D1c_DD.dat'
Rh1c_files = './Results/*Rh1c_DD.dat'

D1c = []
Rh1c = []

#Import files with diffusion constants and Rh values for 1 component model.
for name in glob.glob(path + "*D1c_DD.dat"):
    name = os.path.basename(name)
    D1c = np.append(D1c, name)
    
for name in glob.glob(path + "*Rh1c_DD.dat"):
    name = os.path.basename(name)
    Rh1c = np.append(Rh1c, name)

#Sort files
D1c.sort()
Rh1c.sort()

#Plot D constants
D1 = Import_data(path, D1c[0], D1c[1])
D2 = Import_data(path, D1c[2], D1c[3])
D3 = Import_data(path, D1c[4], D1c[5])
D4 = Import_data(path, D1c[6], D1c[7])        #Duplicate since standard was only measured once.

label = ["frac9", "frac11", "frac12", "sdp_a488"]
ylabel = 'D (um^2/s)'

Plot_bars(D1, D2, D3, D4, label, ylabel, 'Diffusion coefficient', 'diffusion_coefficient', 0)

#Plot Rh
Rh1 = Import_data(path, Rh1c[0], Rh1c[1])
Rh2 = Import_data(path, Rh1c[2], Rh1c[3])
Rh3 = Import_data(path, Rh1c[4], Rh1c[5])
Rh4 = Import_data(path, Rh1c[6], Rh1c[7])    #Duplicate since standard was only measured once.

label = ["frac9", "frac11", "frac12", "sdp_a488"]
ylabel = 'Rh (nm)'

Plot_bars(Rh1, Rh2, Rh3, Rh4, label, ylabel, 'Rh', 'Rh', 1)

##https://www.geeksforgeeks.org/plotting-multiple-bar-charts-using-matplotlib-in-python/?fbclid=IwAR3yOEpVAQz_smHQkzQEfnbuek5AgovowG2db6qcuyhJ1TKVqATCtZkztZE
#how to calculate average error
#move drawing function to separate file

