#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 16:15:41 2021

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
        ax.text(rect.get_x() + rect.get_width()/2., height+(yerr[barCount]/0.95),
                '%.2f' % (height), size=18, color=color, zorder=3,
                ha='center', va='bottom')
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
    w = 0.8    # bar width
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
    
    spread = 0.6    # spread of points as a fraction of bar size
    space = (1 - spread)/2    #extra space on one side of bar
    for i in range(len(x)):
        # distribute scatter randomly across whole width of bar
        ax.scatter(x[i] + np.random.random(y[i].size) * (w*spread) - (w / 2) + w*space, y[i], color="black", s=100, marker='*', zorder=5)

    ax.set_ylim(-0,ylim)
    
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(3)
        
    # increase tick width
    ax.tick_params(width=2)

    return fig


path = './Results/'

''''''''''''''''''''''''''''''''''''
'''  CHANGE INPUT OPTIONS HERE   '''
''''''''''''''''''''''''''''''''''''
# Filenames and calibration error
name = ["hiapp_Na594_488_121721_21-1C", 'rat_iapp_Na594_488_121721_21-1C']
rel_err = 40/350     # error from calibration

#Plot options
plot_title = "Human and Rat IAo_SDS"
strings = np.array(['_D1c_', '_Rh1c_', '_N1c_'])    #Can be expanded if needed.

y_labels = np.array([r'$D$ ($\mathrm{\mu m^2 s^{-1})}$', r'$R_{h}$ (nm)', '# molecules'])
ylim = np.array([350,5,10])    #D, Rh, N ...

bar_labels = np.array(['human_a488_594', 'rat_a488_594'])
save_text = 'human_and_rat_iapp_a594_oligomer'       #text used for file naming

# Donor/Acceptor
key = 'AA'
''''''''''''''''''''''''''''''

##########################
### Read and Plot Data ###
##########################

for i in range(0,len(strings)):
    #Read Data
    # Human IAPP
    t1_h1, h1, err_h1 = rd.ReadData(path + name[0] + strings[i] + key + '.dat')
    #t1_h2, h2, err_h2 = rd.ReadData(path + name[1] + strings[i] + key + '.dat')
    
    # Rat IAPP
    t1_r1, r1, err_r1 = rd.ReadData(path + name[1] + strings[i] + key + '.dat')
    #t1_r2, r2, err_r2 = rd.ReadData(path + name[3] + strings[i] + key + '.dat')
    
    # Plot Data
    x = np.array(range(1, len(name)+1)) # x-coordinates of your bars
    y = [h1, r1]     # data series
    
    filename = save_text + strings[i] + key + '_bar_plot.pdf'
    fig = plotBars(x,y,ylim[i], bar_labels, filename)
    
    # Edit axes, titles, and textsizes
    plt.title(plot_title, fontsize=24)
    plt.ylabel(y_labels[i], fontsize=24)
    plt.rc('font', family='serif',serif='Times')
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=24)
    plt.rc('axes',labelsize=24)
    
    #Save plot
    fig.savefig("./Results/" + filename, dpi = 300, bbox_inches='tight')
    
    

    
