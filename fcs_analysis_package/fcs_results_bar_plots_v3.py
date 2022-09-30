#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 16:15:41 2021

- Update - Jan 17 2022
Will now automatically detect number of columns and adjust plot accordingly.

@author: bab226

Purpose: Secondary analysis of Batch_fitting results of FCS data. Will 
plot bar graphs superimposed with individual data points.

"""

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("./lib") #point python to location containing the below three modules
import readData as rd

# Custom functions
def autolabel(rects, ax):
    """
    Attach a text label above each bar displaying its height.
    
    Parameters
    ----------
    rects : a bar graph
    ax : subplot to add labels

    Returns
    -------
    None.

    """
    barCount = 0 
    colors = plt.cm.cool(np.linspace(0.2,0.7,len(x)))    # Set colormap
    for rect in rects:
        yerr=[((rel_err)*np.mean(yi)+np.std(yi)) for yi in y]
        color = colors[barCount]
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height+(yerr[barCount]/0.90),
                '%.2f' % (height), size=16, color=color, zorder=3,
                ha='center', va='bottom', rotation=75)
        barCount = barCount + 1

def plotBars(x,y,ylim,bar_labels,filename):
    """
    Plots bar graph of data along with scatter plot of individual points.

    Parameters
    ----------
    x : x-coordinates of bars.
    y : data series.
    ylim : option to fix y axis range for consistency between data sets.
    bar_labels : labels for bars.
    filename : text to use as filename.
    
    Returns
    -------
    Figure of Bar Graph.

    """
    
    # Specifiy dimentions of plot, colors, data to include.
    width = 8
    height = width / 1.2
    w = 0.6    # bar width
    colors = plt.cm.cool(np.linspace(0.2,0.7,len(x)))    # Set colormap        

    # Construct bar graph
    fig, ax = plt.subplots()
    
    rect1 = ax.bar(x,
            height=[np.mean(yi) for yi in y],
            yerr=[((rel_err)*np.mean(yi)+np.std(yi)) for yi in y],    # error bars
            capsize=20, # error bar cap width in points
            width=w,    # bar width
            tick_label=bar_labels,
            
            color=(0.5,0.5,0.5,0.3),  # face color transparent
            edgecolor=colors,
            linewidth=4,
            zorder=1
            )
    
    # Change size of plot based on settings above.
    fig.set_size_inches(width, height)
    
    autolabel(rect1, ax)
    
    ax.set_xticklabels(bar_labels, rotation = 75)
    
    spread = 0.6    # spread of points as a fraction of bar size
    space = (1 - spread)/2    #extra space on one side of bar
    for i in range(len(x)):
        # distribute scatter randomly across whole width of bar
        ax.scatter(x[i] + np.random.random(np.size(y[i])) * (w*spread) - (w / 2) + w*space, y[i], color="black", s=10, marker='*', zorder=5)

    ax.set_ylim(-0,ylim)
    
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3)
        
    # increase tick width
    ax.tick_params(width=2)

    return fig



''''''''''''''''''''''''''''''''''''
'''  CHANGE INPUT OPTIONS HERE   '''
''''''''''''''''''''''''''''''''''''
# Filenames and calibration error
name = ['0011x_sds_02-11-22', '0017x_sds_02-11-22', '0026x_sds_02-11-22','0031x_sds_02-11-22',
         '0037x_sds_02-11-22', '0056x_sds_02-11-22', '0083x_sds_02-11-22', '0125x_sds_02-11-22', '025x_sds_02-11-22']
# name = ['0011x_sds_20h_02-12-22', '0017x_sds_20h_02-12-22', '0026x_sds_20h_02-12-22','0031x_sds_20h_02-12-22',
#           '0037x_sds_20h_02-12-22', '0056x_sds_20h_02-12-22', '0083x_sds_20h_02-12-22', '0125x_sds_20h_02-12-22', '025x_sds_20h_02-12-22']
# name = ['sds_iao_dilution_01-28-22_00', 'sds_iao_dilution_01-28-22_01', 'sds_iao_dilution_01-28-22_02', 'sds_iao_dilution_01-28-22_03',
#         'sds_iao_dilution_01-28-22_04', 'sds_iao_dilution_01-28-22_05', 'sds_iao_dilution_01-28-22_06', 'sds_iao_dilution_01-28-22_07',
#         'sds_iao_dilution_01-28-22_08', 'sds_iao_dilution_01-28-22_09']
#name = ['1a_03-01-22', '1b_03-01-22', '2a_03-01-22', '2b_03-01-22', '3a_03-01-22', '3b_03-01-22', '1a_dialysis_12h_03-02-22', '1b_dialysis_12h_03-02-22', 
#        '2a_dialysis_12h_03-02-22', '2b_dialysis_12h_03-02-22', '3a_dialysis_12h_03-02-22', '3b_dialysis_12h_03-02-22']
#name = ['05ul-2_t0h_03-29-22','05ul-2_t2h_03-29-22','05ul-2_t4h_03-29-22', '05ul-2_t12h_03-30-22']

rel_err = 40/470     # error from calibration

#Plot options
# plot_title = "SDS"
# strings = np.array(['D1c_', 'Rh1c_', 'N1c_', 'Rh2c_', 'p1_'])    #Can be expanded if needed.
strings = np.array(['D2c_slow_', 'Rh2c_', 'N2c_', 'p1_']) 
# y_labels = np.array([r'$D$ ($\mathrm{\mu m^2 s^{-1})}$', r'$R_{h}$ (nm)', '# molecules', r'$R_{h}^{slow}$ (nm)', r'$f_{slow}$'])
y_labels = np.array([r'$D$ ($\mathrm{\mu m^2 s^{-1})}$', r'$R_{h}^{slow}$ (nm)', '# molecules', r'$f_{slow}$'])
ylim = np.array([300,6,5,0.8])    #[D, Rh, N ...]

bar_labels = np.array(['t = 0h','t = 2h','t = 4h','t = 12h'])
save_text = 'sds_iao_multipoint_dialysis_50uM_2_'       #text used for file naming

# Donor/Acceptor
key = 'DD'

# MW of SDS
mw = 288.372
''''''''''''''''''''''''''''''

##########################
### Read and Plot Data ###
##########################

''' Bar Graphs '''
for i in range(0,len(strings)):
    #This loop collects data from each sample, for each variable (D, Rh, N, etc.)
    #Array to hold data
    y = []
    
    #Inner loop is for each sample.
    for s in range(0,len(name)):
        # Pool data from each sample
        path = '/Users/bab226/Documents/yale_research/iapp/fcs/fcs_analyzer_package/fcs_analysis_package/Figures/' + name[s] + '/'
        t1, s1, err_s1 = rd.ReadData(path + strings[i] + key + '.dat')
        y.append([s1])
    
    x = np.array(range(0, len(name))) # x-coordinates of your bars
    y = y 
    
    filename = save_text + strings[i] + key + '_bar_plot.pdf'
    fig = plotBars(x,y,ylim[i], bar_labels, filename)
    
    # Edit axes, titles, and textsizes
    # plt.title(plot_title, fontsize=30)
    plt.ylabel(y_labels[i], labelpad=(10), fontsize=24)
    plt.rc('font', family='serif',serif='Times')
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=24)
    plt.rc('axes',labelsize=24)
    
    #Save plot
    path = '/Users/bab226/Documents/yale_research/iapp/fcs/fcs_analyzer_package/fcs_analysis_package/Figures/'
    fig.savefig(path + filename, dpi = 300, bbox_inches='tight')
exit
''' Construct scatter plots '''
''' Rh '''

#This loop collects data from each sample, for each variable (D, Rh, N, etc.)
#Array to hold data
rh2 = []
rh1 = []
#Inner loop is for each sample.
for s in range(0,len(name)):
    # Pool data from each sample
    path = '/Users/bab226/Documents/yale_research/iapp/fcs/results/BB_thermo_sds_iao_step_dilutions_01-28-22/' + name[s] + '/'
    t1, s1, err_s1 = rd.ReadData(path + 'Rh2c_DD' + '.dat')
    rh2.append([s1])
    
for s in range(0,len(name)):
    # Pool data from each sample
    path = '/Users/bab226/Documents/yale_research/iapp/fcs/results/BB_thermo_sds_iao_step_dilutions_01-28-22/' + name[s] + '/'
    t1, s1, err_s1 = rd.ReadData(path + 'Rh1c_DD' + '.dat')
    rh1.append([s1])
    
# Avergage Data
rh1_mean = [np.mean(rh1i) for rh1i in rh1]
rh2_mean = [np.mean(rh2i) for rh2i in rh2]
yerr1=[((rel_err)*np.mean(rh1i)+np.std(rh1i)) for rh1i in rh1]
yerr2=[((rel_err)*np.mean(rh2i)+np.std(rh2i)) for rh2i in rh2]

colors = plt.cm.cool(np.linspace(0.2,0.7,len(x)))    # Set colormap  
fig = plt.figure()
ax = fig.add_subplot()
x = np.array([1.0, 0.5, 0.25, 0.125, 0.083, 0.056, 0.037, 0.031, 0.026, 0.017, 0.011])
ax.errorbar(x[0:len(x)], rh2_mean[0:len(x)], yerr=yerr2[0:len(x)], xerr=0.01, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 5, color='black')
ax.set_xlabel('%SDS (w/v)', labelpad=10)
ax.set_ylim(-0,8)
# Configure second x-axis
x2 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
sds_conc = ((x2/100)/mw)*1000**2  # [SDS] in mM
mic_conc = (sds_conc/62)*1000
mic_conc = np.round(mic_conc, 0)

ax2 = ax.twiny()
ax2.set_xticks(mic_conc)
ax2.set_xlabel(u"SDS micelle conc (\u03bcM)", labelpad=10)
# ax.set_xscale('log')
# change all spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)
        
# increase tick width
ax.tick_params(width=2)
    
ax.set_ylabel(r'$R_{h}^{slow}$ (nm)', labelpad=10)
# plt.savefig('/Users/bab226/Documents/yale_research/iapp/fcs/results/BB_thermo_sds_iao_step_dilutions_01-28-22/' + save_text + 'Rh2c_vs_sds_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')

''' fraction slow '''

# Construct scatter plot
#This loop collects data from each sample, for each variable (D, Rh, N, etc.)
#Array to hold data
y = []
#Inner loop is for each sample.
for s in range(0,len(name)):
    # Pool data from each sample
    path = '/Users/bab226/Documents/yale_research/iapp/fcs/results/BB_thermo_sds_iao_step_dilutions_01-28-22/' + name[s] + '/'
    t1, s1, err_s1 = rd.ReadData(path + 'p1_DD' + '.dat')
    y.append([s1])
    
# Avergage Data
y_mean = [np.mean(yi) for yi in y]
yerr=[((rel_err)*np.mean(yi)+np.std(yi)) for yi in y]

colors = plt.cm.cool(np.linspace(0.2,0.7,len(x)))    # Set colormap  
fig, ax = plt.subplots()
x = np.array([1.0, 0.5, 0.25, 0.125, 0.083, 0.056, 0.037, 0.031, 0.026, 0.017])
ax.errorbar(x[0:len(x)], y_mean[0:len(x)], yerr=yerr[0:len(x)], xerr=0.01, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 5, color='black')
ax.set_xlabel('%SDS (w/v)', labelpad=10)
ax.set_ylim(-0,1)
# Configure second x-axis
x2 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
sds_conc = ((x2/100)/mw)*1000**2  # [SDS] in mM
mic_conc = (sds_conc/62)*1000
mic_conc = np.round(mic_conc, 0)
ax2 = ax.twiny()
ax2.set_xticks(mic_conc)
ax2.set_xlabel(u"SDS micelle conc (\u03bcM)", labelpad=10)
# ax.set_xscale('log')
# change all spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)
        
# increase tick width
ax.tick_params(width=2)

ax.set_ylabel(r'$f_{slow}$', labelpad=10)
# plt.savefig('/Users/bab226/Documents/yale_research/iapp/fcs/results/BB_thermo_sds_iao_step_dilutions_01-28-22/' + save_text + 'p1_vs_sds_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')

''' Concentration '''

# Construct scatter plot
#This loop collects data from each sample, for each variable (D, Rh, N, etc.)
#Array to hold data
y = []
#Inner loop is for each sample.
for s in range(0,len(name)):
    # Pool data from each sample
    path = '/Users/bab226/Documents/yale_research/iapp/fcs/results/BB_thermo_sds_iao_step_dilutions_01-28-22/' + name[s] + '/'
    t1, s1, err_s1 = rd.ReadData(path + 'C2c_DD' + '.dat')
    y.append([s1])
    
# Avergage Data
y_mean = [np.mean(yi) for yi in y]
yerr=[((rel_err)*np.mean(yi)+np.std(yi)) for yi in y]

colors = plt.cm.cool(np.linspace(0.2,0.7,len(x)))    # Set colormap  
fig, ax = plt.subplots()
x = np.array([1.0, 0.5, 0.25, 0.125, 0.083, 0.056, 0.037, 0.031, 0.026, 0.017])
ax.errorbar(x[0:len(x)], y_mean[0:len(x)], yerr=yerr[0:len(x)], xerr=0.01, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 5, color='black')
ax.set_xlabel('%SDS (w/v)', labelpad=10)
ax.set_ylim(-0,40)
# Configure second x-axis
x2 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
sds_conc = ((x2/100)/mw)*1000**2  # [SDS] in mM
mic_conc = (sds_conc/62)*1000
mic_conc = np.round(mic_conc, 0)
ax2 = ax.twiny()
ax2.set_xticks(mic_conc)
ax2.set_xlabel(u"SDS micelle conc (\u03bcM)", labelpad=10)
# ax.set_xscale('log')
# change all spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(3)
        
# increase tick width
ax.tick_params(width=2)
ax.set_ylabel('Concentration (nM)', labelpad=10)
# plt.savefig('/Users/bab226/Documents/yale_research/iapp/fcs/results/BB_thermo_sds_iao_step_dilutions_01-28-22/' + save_text + 'c2_vs_sds_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')





















