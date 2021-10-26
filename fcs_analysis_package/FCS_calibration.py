# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 16:13:06 2021

@author: gwg24 (Greg Gomes)
"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import lmfit
from uncertainties import ufloat
from uncertainties.umath import *

import sys
sys.path.append("./lib") #point python to location containing the below three modules
import FCS_fitfunc as ff
import SPT_reader_edit as spt
import FCS_helpful as fcs


# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1


''' Functions to facilitate  fitting
'''
def calibration(t, G, err, kappa = 5):
    # 1-comp diffusion
    
    A0_guess = np.mean(G[:5])
    arg = abs(G - A0_guess/2)
    ind = np.argmin(arg)
    tau_guess = t[ind]    
    print('Guessing tau = %f ms' %(tau_guess))
    
    model = lmfit.Model(ff.diffusion_3d)
    params = model.make_params(A0=A0_guess, tau_diff=tau_guess)
    params['A0'].set(min=0.00001, value=A0_guess)
    params['tau_diff'].set(min=1e-6, value=tau_guess)
    params['Ginf'].set(value=0, vary = False)
    params['kappa'].set(value=kappa, vary=True)  # 3D model only
    
    weights = 1/err
    t_max = 1e7
    fitres = model.fit(G[t < t_max], timelag=t[t<t_max], params=params, method='least_squares',
                       weights=weights[t<t_max])
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
      
    return fitres



# measurement_group = spt.Read_FCS('./Data/A488_A594_FRET_FCS_grouped.dat')

# measurement_group = spt.Read_FCS('/Users/bab226/Documents/yale_research/iapp/fcs/fcs-analysis-package/Data/BB_IAPP_SEC_fractions_08-17-21.sptw/calibration/')
measurement_group = spt.Read_FCS('/Users/bab226/Documents/yale_research/iapp/fcs/fcs_analyzer_package/fcs_analysis_package/Data/BB_dextran_mixtures.sptw/calibration/rho110_10_20_21')

# Possible key values are DD (autocorrelation Donor Donor), AA (auto, accepptor acceptor), DxA (donor acceptor cross correlation)
key = 'DD'

#Intialize lists to store results of fitting, 2c  = 2 component, 1c = 1 component.

tau_1c = []
tau_err1c = []
N_1c = []
N_1c_err = []
kappa = []
kappa_err = []
t_measurement = [] #measurement time which will be grepped from title of measurement
#If plotting time series on same graph, this can be helpful
#There are n measurements in measurement group, so we divide a color map into n linearly spaced colors
#i.e., plot goes from hot to cold colors or whatever --can change color map from jet to other
ncol = len(measurement_group)
colors = plt.cm.jet(np.linspace(0,1,ncol))

i = 0 

for m in measurement_group:
    #SPT64 saves file name in measurement group as e.g., 'ThisName_T1800s_1.ptu'
    #where T1800s means this measurement was taken at 1800 seconds after start.
    #Use regular expression to find pattern Tsomeinteger and then convert it to a number
    try:
        measurement_time = float(re.findall(r'T\d+',m['name'])[0][1:]) 
    except:
        measurement_time = 0    
    measurement_time = (measurement_time/60) #now minutes 
    mylabel = 't = %.2f min' %(measurement_time)

    #Use key to decide which correlation curve to look at
    t = m[key]['time']
    G = m[key]['G']
    err = m[key]['err']
 

    print()
    print()
    print('############')
    print('Results of measurement at time t = %.2f min' %(measurement_time))
    fitres1c = calibration(t, G, err)

    
    #Store results
    t_measurement.append(measurement_time)    
    tau_1c.append(fitres1c.values['tau_diff'])
    tau_err1c.append(fitres1c.params['tau_diff'].stderr)
    A01c = ufloat(fitres1c.values['A0'], fitres1c.params['A0'].stderr) 
    x = 1/A01c #N = 1/A
    N_1c.append(x.nominal_value)
    N_1c_err.append(x.std_dev)
    kappa.append(fitres1c.values['kappa'])
    kappa_err.append(fitres1c.params['kappa'].stderr)
    
    #Plot each separately
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t,G, yerr = err, linewidth =1, label = mylabel, linestyle = '', marker = 'o', markersize = 2, color = colors[i])
    
    
    ax.plot(t, fitres1c.best_fit, color = 'k', label = '1-comp', linestyle = '-', zorder = 10) #zorder forces best fit to plot on top errorbar
    
    ax.set_xscale('log')
    ax.set_xlabel(r'$\tau$ (ms)', labelpad=10)
    ax.set_ylabel(r'$G(\tau)$', labelpad=10)
    ax.legend()
    # savefig_string = './Figures/' + m['name'][:-4] + '_' + key + '.png' #saves in Figures folder, as measurement name, appended with which correlation
    # plt.savefig(savefig_string, dpi=300, transparent=False, bbox_inches='tight')
    i += 1
    
    
tau_1c = np.array(tau_1c)
tau_err1c = np.array(tau_err1c)
N_1c = np.array(N_1c)
N_1c_err = np.array(N_1c_err)
N_1c = np.array(N_1c)
N_1c_err = np.array(N_1c_err)
kappa = np.array(kappa)


print()
print()
print('############')
print('Quick summary')
print('Analyzing ' + key + ' curves')
print('Mean 1-comp diffusion time = %.4f +/- %.5f ms' %(np.mean(tau_1c), np.std(tau_1c, ddof =1)))
print('Mean 1-comp kappa = %.4f +/- %.5f ' %(np.mean(kappa), np.std(kappa, ddof =1)))


#color the plots according to whether you are looking at DD,AA or DxA
color_dict = {"DD": 'b',
              "AA": 'r',
              "DxA": 'k'
              }


''' Plots for 1-comp fitting 
'''''''''''''''''''''''''''''''''''''''
'''
'''
width = 3.42
fig = plt.figure(figsize=(width,width/1.62))
ax = fig.add_axes([0, 0, 1, 1])
ax.errorbar(t_measurement, tau_1c*1000, yerr = tau_err1c*1000, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
ax.set_xlabel(r'time (min)', labelpad=10)
ax.set_ylabel(r'$t_{d}$ ($\mathrm{\mu s}$)', labelpad=10)
plt.savefig('./Figures/comp1_diffusion_time_' +key + '.png', dpi=300, transparent=False, bbox_inches='tight')

width = 3.42
fig = plt.figure(figsize=(width,width/1.62))
ax = fig.add_axes([0, 0, 1, 1])
ax.errorbar(t_measurement, N_1c, yerr = N_1c_err, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
ax.set_xlabel(r'time (min)', labelpad=10)
ax.set_ylabel(r'N', labelpad=10)
plt.savefig('./Figures/Nmol1comp_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')


width = 3.42
fig = plt.figure(figsize=(width,width/1.62))
ax = fig.add_axes([0, 0, 1, 1])
ax.errorbar(t_measurement, kappa, yerr = kappa_err, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
ax.set_xlabel(r'time (min)', labelpad=10)
ax.set_ylabel(r'kappa', labelpad=10)
plt.savefig('./Figures/kappa_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')

