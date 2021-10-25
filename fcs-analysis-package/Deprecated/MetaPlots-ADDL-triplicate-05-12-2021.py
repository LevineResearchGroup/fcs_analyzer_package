# -*- coding: utf-8 -*-
"""
Created on Wed May 12 11:48:30 2021

@author: gwg24
"""

import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# from uncertainties import ufloat
# from uncertainties.umath import *

from scipy import stats
import seaborn as sns

mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 15
plt.rcParams['axes.linewidth'] = 2
plt.rcParams["errorbar.capsize"] = 4
lw = 2

colors = sns.color_palette("Set2", 3 )
# colors = sns.color_palette("rocket", 3 )
color1 = colors[0]
color2 = colors[1]
color3 = colors[2]

color1 = (27/255,158/255,119/255)
color2 = (217/255, 95/255, 2/255)
color3 = (117/255, 112/255, 179/255)

# color2 = 'palevioletred'
# color1 = 'lightskyblue'

#%%
def ReadData(filename):
    
    '''
    WARNING: I ADDED A LIMITER (ONLY RETURN FIRST FIVE) BECAUSE I WANTED EACH TO HAVE ONLY 5 DATA POINTS!
    '''
    sample_size = 3
    t = []
    R = []
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
            R.append(float(row[1]))
            err.append(float(row[2]))


    return np.array(t)[:sample_size], np.array(R)[:sample_size], np.array(err)[:sample_size]


def PlotTubeGrouop(var1, var2, ylabel, plottitle, savetitle):
    
    #var1,2,3 should be Tube PBS, Tube SDS, and Tube X in this case (see labels)
    n = 3 #number of samples --> for pltting std of the mean ( std/ srt(n)) rather than SD
    
    labels = ['A', 'B']
    width = 0.25  # the width of the bars
    # Set position of bar on X axis
    r1 = np.arange(len(labels))
    r2 = [x + width for x in r1]
    # r3 = [x - width for x in r1]
    figwidth = 3.42*2
    fig = plt.figure(figsize=(figwidth,figwidth/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    rects1 = ax.bar(r1, var1[:,0], yerr = var1[:,1]/np.sqrt(n), width = 0.8*width, label='PBS', color = color1, alpha = 1, edgecolor = 'k', linewidth = lw)
    rects2 = ax.bar(r2, var2[:,0], yerr = var2[:,1]/np.sqrt(n), width = 0.8*width, label='SDS', color = color2, alpha = 1, edgecolor = 'k', linewidth = lw)
    # rects3 = ax.bar(r3, var3[:,0], yerr = var3[:,1]/np.sqrt(n), width = 0.8*width, label='TX', color = color3, alpha = 1, edgecolor = 'k', linewidth = lw)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Sample')
    ax.set_title(plottitle)
    ax.set_xticks(r1)
    ax.set_xticklabels(labels)
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    ax.bar_label(rects1, padding=2, fmt = '%.2f')
    ax.bar_label(rects2, padding=2, fmt = '%.2f')
    # ax.bar_label(rects3, padding=3, fmt = '%.2f')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    
    fig.tight_layout()
    plt.savefig(savetitle, dpi=300, transparent=False, bbox_inches='tight')
    
    plt.show()


def ttest_custom(x,y,n):
    
    x_means = x[:,0]
    x_std = x[:,1]
    y_means = y[:,0]
    y_std = y[:,1]
    sDelta = np.sqrt( (x_std**2)/n + (y_std**2)/n  ) #Welch's t test, unequal variances, unequal sample sizes
    t = (x_means - y_means)/sDelta
    df = ((sDelta**2)**2) / (  (((x_std**2)/n)**2)/(n-1) + (((y_std**2)/n)**2)/(n-1)     )
    
    for tt,ddf in zip(t,df):
        pval = stats.t.sf(np.abs(tt), ddf)*2  # two-sided pvalue = Prob(abs(t)>tt)
        print('t-statistic = %6.3f pvalue = %6.4f' % (tt, pval))

def ttest_byaliquot(x,n):
#
    x_means = x[:,0]
    x_std = x[:,1]
    for i in range(len(x_means)):
        for j in range(len(x_means)):
            sDelta = np.sqrt( (x_std[i]**2)/n + (x_std[j]**2)/n  ) #Welch's t test, unequal variances, unequal sample sizes
            t = (x_means[i] - x_means[j])/sDelta
            df = ((sDelta**2)**2) / (  (((x_std[i]**2)/n)**2)/(n-1) + (((x_std[j]**2)/n)**2)/(n-1)     )
            pval = stats.t.sf(np.abs(t), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
            print('A %d vs A %d t-statistic = %6.3f pvalue = %6.4f' % (i+1,j+1,t, pval))  
    


''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''' 'Load Data' ''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''

#Diffusion
D_frac9 = []
t, R, err = ReadData('./Results/frac9_sample1_D1c_DD.dat')
D_frac9.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/frac9_sample2_D1c_DD.dat')
D_frac9.append([np.mean(R), np.std(R,ddof = 1)])

'''

'''
D_frac11 = []
t, R, err = ReadData('./Results/frac11_sample1_D1c_DD.dat')
D_frac11.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/frac11_sample2_D1c_DD.dat')
D_frac11.append([np.mean(R), np.std(R,ddof = 1)])
# t, R, err = ReadData('./Results/ADDL_Tube82_Aliquot2_Rh2c_DD.dat')
# Rh_T82.append([np.mean(R), np.std(R,ddof = 1)])
# t, R, err = ReadData('./Results/ADDL_Tube82_aliquot3_Rh2c_DD.dat')
# Rh_T82.append([np.mean(R), np.std(R,ddof = 1)])
# t, R, err = ReadData('./Results/ADDL_Tube82_Aliquot4_Rh2c_DD.dat')
# Rh_T82.append([np.mean(R), np.std(R,ddof = 1)])
#%%

D_frac9 = np.array(D_frac9)
D_SDS = np.array(D_frac11)


PlotTubeGrouop(D_frac9, D_frac11, ylabel = 'D (um^s/s)', plottitle = 'Diffusion coef by sample', savetitle = 'D_plot.png')


ttest_custom(D_frac9, D_SDS, n = 3)

# print('Tube 81, Rh Aliquot pair differences')
# ttest_byaliquot(Rh_T81, n=5)
# print('Tube 82, Rh Aliquot pair differences')
# ttest_byaliquot(Rh_T82, n=5)
# print('Tube X, Rh Aliquot pair differences')
# ttest_byaliquot(Rh_TX, n=5)







''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''' 'fraction slow' ''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''
"""
'''
Tube 81 Rh
'''
p_T81 = []
t, R, err = ReadData('./Results/ADDL_Tube81_Aliquot1_group_p1_2c_DD.dat')
p_T81.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_Tube81_Aliquot2_p1_2c_DD.dat')
p_T81.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_Tube81_aliquot3_p1_2c_DD.dat')
p_T81.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_Tube81_Aliquot4_p1_2c_DD.dat')
p_T81.append([np.mean(R), np.std(R,ddof = 1)])

'''
Tube 82 Rh
'''
p_T82 = []
t, R, err = ReadData('./Results/ADDL_Tube82_Aliquot1_group_p1_2c_DD.dat')
p_T82.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_Tube82_Aliquot2_p1_2c_DD.dat')
p_T82.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_Tube82_aliquot3_p1_2c_DD.dat')
p_T82.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_Tube82_Aliquot4_p1_2c_DD.dat')
p_T82.append([np.mean(R), np.std(R,ddof = 1)])

'''
Tube X Rh
'''
p_TX = []
t, R, err = ReadData('./Results/ADDL_TubeUnnumb_Aliquot1_group_p1_2c_DD.dat')
p_TX.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_TubeUnnumb_Aliquot2_p1_2c_DD.dat')
p_TX.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_TubeX_aliquot3_p1_2c_DD.dat')
p_TX.append([np.mean(R), np.std(R,ddof = 1)])
t, R, err = ReadData('./Results/ADDL_TubeX_Aliquot4_p1_2c_DD.dat')
p_TX.append([np.mean(R), np.std(R,ddof = 1)])

p_T81 = np.array(p_T81)
p_T82 = np.array(p_T82)
p_TX = np.array(p_TX)

PlotTubeGrouop(p_T81, p_T82, p_TX, ylabel = r'$f_{slow}$', plottitle = r'$f_{slow}$ by Tube and Aliquot', savetitle = 'fslow_bar_v2.png')


print('# fslow tube differences, by aliquot')
print('#Tube 81 vs 82 ')
ttest_custom(p_T81, p_T82, n = 5)
print('#Tube 81 vs X ')
ttest_custom(p_T81, p_TX, n = 5)
print('#Tube 82 vs X ')
ttest_custom(p_T82, p_TX, n = 5)


print('Tube 81, fslow Aliquot pair differences')
ttest_byaliquot(p_T81, n=5)
print('Tube 82, fslow Aliquot pair differences')
ttest_byaliquot(p_T82, n=5)
print('Tube X, fslow Aliquot pair differences')
ttest_byaliquot(p_TX, n=5)
"""