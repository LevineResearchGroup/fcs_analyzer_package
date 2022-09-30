#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 15:09:59 2022

Wrapper for Rh Data from Greg's script
Export Rh data for csv for easy imput into excell or Prism

@author: bab226
"""

import numpy as np
import gromacs.formats as gmx
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import readData as rd
from uncertainties import ufloat

def get_cmap(n, name='viridis'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

color = sns.color_palette('deep', 5, desat=1)

#Data Path:
path = "/Users/bab226/Documents/yale_research/iapp/fcs/fcs_analyzer_package/fcs_analysis_package/Figures/sds_iao_dilution_02-11-22/"

f_arr = []
err_arr = []
for i in range(0,10):
    i = str(i)
    i = i.zfill(2)  # padding 0 
    directory = 'sds_iao_dilution_02-11-22_%s/' %(i)

    file = 'p1_DD.dat'

    t1, fslow, err = rd.ReadData(path + directory + file, )
    
    mean_f = np.mean(fslow)
    sd_f = np.std(fslow)
    
    f_arr.append(mean_f)
    err_arr.append(sd_f)


df = pd.DataFrame({'fslow': f_arr, 'err': list(err_arr)}, columns=['fslow', 'err'])

df.to_csv(path + 'f.csv', index=False)