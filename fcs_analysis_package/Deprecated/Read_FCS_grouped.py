# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 11:33:46 2021

@author: gwg24
"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import lmfit
from uncertainties import ufloat
from uncertainties.umath import *


import FCS_fitfunc as ff
import SPT_reader as spt
import FCS_helpful as fcs


# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1


'''
1/G0 plots
'''



'''
IN SPT for FRET-FCS
A = Ch2, Al488, 0-50 ns
B = Ch1, Al594, 0-50 ns
Do not export with fits.
'''


''' Functions to facilitate  fitting'''
    
def simplefit(t, G, err):
    # 1-comp diffusion
    
    A0_guess = np.mean(G[:5])
    arg = abs(G - A0_guess/2)
    ind = np.argmin(arg)
    tau_guess = t[ind]    
    print('Guessing tau = %f ms' %(tau_guess))
    
    model = lmfit.Model(ff.diffusion_3d)
    params = model.make_params(A0=A0_guess, tau_diff=tau_guess)
    params['A0'].set(min=0.01, value=A0_guess)
    params['tau_diff'].set(min=1e-6, value=tau_guess)
    params['Ginf'].set(value=0, vary = False)
    params['kappa'].set(value=5.9, vary=False)  # 3D model only
    
    weights = 1/err
    t_max = 1e7
    fitres = model.fit(G[t < t_max], timelag=t[t<t_max], params=params, method='least_squares',
                       weights=weights[t<t_max])
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
      
    return fitres  

def simplefit2(t, G, err):
    # two component diffusion
    model = lmfit.Model(ff.twocomp_diff3d)
    
    A0_guess = np.mean(G[:5])
    arg = abs(G - A0_guess/2)
    ind = np.argmin(arg)
    tau_guess = t[ind]    
    print('Guessing tau = %f ms' %(tau_guess))
    
    params = model.make_params(A0=A0_guess)
    params['A0'].set(min=0.01, value = A0_guess)
    params['tau_diff1'].set(min=1e-6, value = tau_guess) #slow component
    params['tau_diff2'].set(min=1e-6, value = 0.036, vary = False) #fast component, usually fixed
    params['p1'].set(min = 0, max = 1, value = 0.5) #fraction slow
    
    params['Ginf'].set(value=0, vary = False)
    params['kappa'].set(value=5.9, vary=False)  # 3D model only
    
    # weights = np.ones_like(avg_G)
    weights = 1/err
    t_max = 1e7
    fitres = model.fit(G[t < t_max], timelag=t[t<t_max], params=params, method='least_squares',
                       weights=weights[t<t_max])
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
    
    return fitres
    

def FRETfit2(t, G, err):
    # two component FRET-FCS
    model = lmfit.Model(ff.Gdd_model2)
    params = model.make_params(A0=1.6, tau_diff=0.03)
    params['A0'].set(min=0.0, value=0.05)
    params['p1'].set(min = 0, max = 1, value = 0.5)
    params['tau_diff1'].set(min=1e-6, value=0.036, vary = False)
    params['tau_diff2'].set(min=0, value=0.452)
    params['Ginf'].set(value=0, vary = True)
    params['kappa'].set(value=5.9, vary=False)  # 3D model only
    
    params['tau_k'].set(min = 0.0, value=0.036)
    params['a'].set( value = 0.1, min = 0, max = 1 )
    params['b'].set( value = 0, vary =False)    
    # weights = np.ones_like(avg_G)
    weights = 1/err
    t_max = 1e7
    fitres = model.fit(G[t < t_max], timelag=t[t<t_max], params=params, method='least_squares',
                       weights=weights[t<t_max])
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
    
    # print(fitres.params['tau_diff'].stderr * 1e6)
    
    if fitres.params['tau_diff2'].stderr is not None:
        stderr_us = fitres.params['tau_diff2'].stderr
    else:
        stderr_us = fitres.values['tau_diff2']*100
    
    return fitres, fitres.values['tau_diff2'], stderr_us 

def FRETfit_global(t, Gdd, errdd, Gaa, erraa, Gx, errx):
        
    global_data, global_err = ff.global_E_data(Gdd,Gaa,Gx, errdd, erraa, errx)
    
    model = lmfit.Model(ff.global_E_model)
    params = model.make_params(A0=1.6, tau_diff=100e-6)
    params['A0'].set(min=0.0, value=0.05)
    params['tau_diff'].set(min=0, value=0.452)
    params['Ginf'].set(value=0, vary = False)
    params['kappa'].set(value=5.9, vary=False)  # 3D model only
    
    params['tau_k'].set(min = 0.0, value=1.6e-3)
    params['a'].set(min = 0.0, value = 0.1 )
    params['b'].set(min = 0.0, value = 0.1)
    params['c'].set(min = 0.0, value = 0.1 )
    params['d'].set(min = 0.0, value = 0.1)      
    # weights = np.ones_like(avg_G)
    weights = 1/global_err
    fitres = model.fit(global_data, timelag=t, params=params, method='least_squares',
                       weights=weights)
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
    return fitres

def FRETfit_global2(t, Gdd, errdd, Gaa, erraa, Gx, errx):
        
    global_data, global_err = ff.global_E_data(Gdd,Gaa,Gx, errdd, erraa, errx)
    
    model = lmfit.Model(ff.global_E_model2)
    params = model.make_params(A0=1.6)
    params['A0'].set(min=0.0, value=0.5)
    params['p1'].set(min = 0, max = 1, value = 0.5)
    params['tau_diff1'].set(min=0, value=0.036, vary = False)
    params['tau_diff2'].set(min=0, value=0.452)
    params['Ginf'].set(value=0, vary = False)
    params['kappa'].set(value=5.9, vary=False)  # 3D model only
    
    params['tau_k'].set(min = 0.0, value=1.6e-3)
    params['a'].set(min = 0.0, value = 0.1 )
    params['b'].set(min = 0.0, value = 0.1)
    params['c'].set(min = 0.0, value = 0.1 )
    params['d'].set(min = 0.0, value = 0.1)      
    # weights = np.ones_like(avg_G)
    weights = 1/global_err
    fitres = model.fit(global_data, timelag=t, params=params, method='least_squares',
                       weights=weights)
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
    return fitres


# measurement_group = spt.Read_FCS('./Data/A488_A594_FRET_FCS_grouped.dat')
# measurement_group = spt.Read_FCS('./Data/monomericAB_FRET_FCS_grouped.dat')
# measurement_group = spt.Read_FCS('./Data/Barghorn_4xdil_nextday_FRET-FCS.dat')
# measurement_group = spt.Read_FCS('./Data/Barghorn_4xdil_nextday_filt1_FRET-FCS.dat')
# measurement_group = spt.Read_FCS('./Data/Barghorn_4xdil_nextday_4xfilt_FRET-FCS.dat')
# measurement_group = spt.Read_FCS('./Data/Barghorn_4xdil_FRET-FCS.dat')
# measurement_group = spt.Read_FCS('./Data/Barghorn_weekend_FRET-FCS.dat')
measurement_group = spt.Read_FCS('./Data/SDS_excess_ABO-FRET-FCS.dat')
# measurement_group = spt.Read_FCS('')
# measurement_group = spt.Read_FCS('')


if 0:
    
    j = 3
    t = measurement_group[j]['DD']['time']    
    Gdd = measurement_group[j]['DD']['G']    
    errdd = measurement_group[j]['DD']['err']
    taa = measurement_group[j]['AA']['time']    
    Gaa = measurement_group[j]['AA']['G']    
    erraa = measurement_group[j]['AA']['err']
    tx = measurement_group[j]['DxA']['time']
    Gx = measurement_group[j]['DxA']['G']    
    errx = measurement_group[j]['DxA']['err']
    
    # errx = 10000*np.ones(np.size(errx))
    
    # print(len(t))
    # print(len(taa))
    # print(len(tx))
    # fitres_fret = FRETfit_global(t, Gdd, errdd, Gaa, erraa, Gx, errx)
    fitres_fret = FRETfit_global2(t, Gdd, errdd, Gaa, erraa, Gx, errx)

    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t,Gdd, yerr = errdd, linewidth =2, label = 'data', linestyle = '', marker = 'o', markersize = 1, color = 'k')
    ax.plot(t, fitres_fret.best_fit[0:len(t)], color = 'b', label = '3d-diffusion + kinetics')
    ax.set_xscale('log')
    plt.title('Gdd')
    ax.legend()
    
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t,Gaa, yerr = erraa, linewidth =2, label = 'data', linestyle = '', marker = 'o', markersize = 1, color = 'k')
    ax.plot(t, fitres_fret.best_fit[len(t):2*len(t)], color = 'r', label = '3d-diffusion + kinetics')
    plt.title('Gaa')
    ax.set_xscale('log')
    ax.legend()
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t,Gx, yerr = errx, linewidth =2, label = 'data', linestyle = '', marker = 'o', markersize = 1, color = 'k')
    ax.plot(t, fitres_fret.best_fit[2*len(t):3*len(t)], color = 'm', label = '3d-diffusion + kinetics')
    plt.title('Gda')
    ax.set_xscale('log')
    ax.legend()

if 1:
    '''Example of batch fitting/plotting'''
    
    
    # Possible key values are DD (autocorrelation Donor Donor), AA (auto, accepptor acceptor), DxA (donor acceptor cross correlation)
    key = 'DxA'
    
    Toffset = 0 #offset the time axis by some amount (make sure consistent with units for measurement_time (i.e., minutes, hours?))

    
    #Intialize lists to store results of fitting, 2c  = 2 component, 1c = 1 component.
    tau_slow2c = [] 
    tau_slow_err2c = []
    p1 = []
    p1_err = []
    
    tau_1c = []
    tau_err1c = []
    
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
        measurement_time = float(re.findall(r'T\d+',m['name'])[0][1:]) 
        measurement_time = Toffset + (measurement_time/60) #now minutes 
        mylabel = 't = %.2f min' %(measurement_time)

        #Use key to decide which correlation curve to look at
        t = m[key]['time']
        G = m[key]['G']
        err = m[key]['err']
 
               
        print()
        print()
        print('############')
        print('Results of measurement at time t = %.2f min' %(measurement_time))
        fitres1c = simplefit(t, G, err)
        fitres2c = simplefit2(t, G, err)
        
        #Store results
        t_measurement.append(measurement_time)
        tau_slow2c.append(fitres2c.values['tau_diff1'])
        tau_slow_err2c.append(fitres2c.params['tau_diff1'].stderr)
        p1.append(fitres2c.values['p1'])
        p1_err.append(fitres2c.params['p1'].stderr)
        
        tau_1c.append(fitres1c.values['tau_diff'])
        tau_err1c.append(fitres1c.params['tau_diff'].stderr)
        
        
        #Plot each separately
        width = 3.42
        fig = plt.figure(figsize=(width,width/1.62))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.errorbar(t,G, yerr = err, linewidth =1, label = mylabel, linestyle = '', marker = 'o', markersize = 2, color = colors[i])
        
        
        ax.plot(t, fitres2c.best_fit, color = 'k', label = '2-comp', zorder = 10) #zorder forces best fit to plot on top errorbar
        ax.plot(t, fitres1c.best_fit, color = 'm', label = '1-comp', linestyle = '--', zorder = 10) #zorder forces best fit to plot on top errorbar
        
        ax.set_xscale('log')
        ax.set_xlabel(r'$\tau$ (ms)', labelpad=10)
        ax.set_ylabel(r'$G(\tau)$', labelpad=10)
        ax.legend()
        savefig_string = './Figures/' + m['name'][:-4] + '_' + key + '.png' #saves in Figures folder, as measurement name, appended with which correlation
        plt.savefig(savefig_string, dpi=300, transparent=False, bbox_inches='tight')
        i += 1
        
    #turn results into np arrays so we can do maths with them
    t_measurement = np.array(t_measurement)
    tau_slow2c = np.array(tau_slow2c)
    tau_slow_err2c = np.array(tau_slow_err2c)
    p1 = np.array(p1)
    p1_err = np.array(p1_err)
    
    tau_1c = np.array(tau_1c)
    tau_err1c = np.array(tau_err1c)
    
    #convert diffusion times to D's and Rh's -- accepts ufloat type for td_ref and D_ref
    D2c, eD2c = fcs.td2D(tau_slow2c, tau_slow_err2c, temperature_lab = 22, td_ref = ufloat(0.0293, 0.002), D_ref = ufloat(400, 0.05*400), temperature_ref = 22 )
    Rh2c, err_Rh2c = fcs.D2Rh(D2c, eD2c, temperature_lab = 22)
    
    D1c, eD1c = fcs.td2D(tau_1c, tau_err1c, temperature_lab = 22, td_ref = ufloat(0.0293, 0.00033), D_ref = ufloat(400, 0.05*400), temperature_ref = 22 )
    Rh1c, err_Rh1c = fcs.D2Rh(D1c, eD1c, temperature_lab = 22)
    
    
    print()
    print()
    print('############')
    print('Quick summary')
    print('Analyzing ' + key + ' curves')
    print('Mean 1-comp diffusion time = %.4f +/- %.5f ms' %(np.mean(tau_1c), np.std(tau_1c, ddof =1)))
    print('Mean 1-comp diffusion coeff slow = %.4f +/- %.5f um^2 s^-1' %(np.mean(D1c), np.std(D1c, ddof =1)))
    print('Mean 1-comp Rh slow = %.4f +/- %.5f nm' %(np.mean(Rh1c), np.std(Rh1c, ddof =1)))
    
    print('Mean 2-comp slow diffusion time = %.4f +/- %.5f ms' %(np.mean(tau_slow2c), np.std(tau_slow2c, ddof =1)))
    print('Mean 2-comp diffusion coeff slow = %.4f +/- %.5f um^2 s^-1' %(np.mean(D2c), np.std(D2c, ddof =1)))
    print('Mean 2-comp Rh slow = %.4f +/- %.5f nm' %(np.mean(Rh2c), np.std(Rh2c, ddof =1)))

    #color the plots according to whether you are looking at DD,AA or DxA
    color_dict = {"DD": 'b',
                  "AA": 'r',
                  "DxA": 'k'
                  }

    ''' Plots for two-comp fitting '''
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t_measurement, tau_slow2c*1000, yerr = tau_slow_err2c*1000, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'$t_{d}$ ($\mathrm{\mu s}$)', labelpad=10)
    plt.savefig('./Figures/Slow_diffusion_time_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t_measurement, D2c, yerr = eD2c, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'$D_{slow}$ ($\mathrm{\mu m^2 s^{-1}}$)', labelpad=10)
    plt.savefig('./Figures/Slow_diffusion_coeff_' + key +'.png', dpi=300, transparent=False, bbox_inches='tight')
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t_measurement, Rh2c, yerr = err_Rh2c, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'$R_{h}^{slow}$ (nm)', labelpad=10)
    plt.savefig('./Figures/Slow_Rh_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t_measurement, p1, yerr = p1_err, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'$f_{slow}$', labelpad=10)
    plt.savefig('./Figures/Slow_frac_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    
    ''' Plots for 1-comp fitting '''
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
    ax.errorbar(t_measurement, D1c, yerr = eD1c, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'$D$ ($\mathrm{\mu m^2 s^{-1}}$)', labelpad=10)
    plt.savefig('./Figures/comp1_diffusion_coeff_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t_measurement, Rh1c, yerr = err_Rh1c, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'$R_{h}$ (nm)', labelpad=10)
    plt.savefig('./Figures/comp1_Rh_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    

if 0:
    
    j = 6
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])

    t = measurement_group[j]['DxA']['time']
    Gx = measurement_group[j]['DxA']['G']
    errx = measurement_group[j]['DxA']['err']
    
    t = measurement_group[j]['AA']['time']
    Gaa = measurement_group[j]['AA']['G']
    erraa = measurement_group[j]['AA']['err']
    
    t = measurement_group[j]['AA']['time']
    Gdd = measurement_group[j]['DD']['G']
    errdd = measurement_group[j]['DD']['err']
    
    rat_dd_da = Gdd/Gx
    
    # ax.errorbar(t,Gdd, yerr = errdd, linewidth =2, label = 'DD', linestyle = '', marker = 'o', markersize = 1, color = 'b')
    # ax.errorbar(t,Gaa, yerr = erraa, linewidth =2, label = 'AA', linestyle = '', marker = 'o', markersize = 1, color = 'r')    
    # ax.errorbar(t,Gx, yerr = errx, linewidth =2, label = 'DxA', linestyle = '', marker = 'o', markersize = 1, color = 'k')
    
    ax.plot(t,rat_dd_da, linestyle = '', marker = 'o', markersize = 1, color = 'k')        
            
    ax.set_xscale('log')
    ax.legend()
        
    ax.set_xlabel(r'$\tau$ (ms)', labelpad=10)
    ax.set_ylabel(r'$G(\tau)$', labelpad=10)
    # ax.set_xlim(1e-3, 10)
    # ax.set_ylim(0,4)
    # ax.set_ylim(1,3)      
    ymax = 1.5*np.mean(rat_dd_da[:5])
    ymin = 0.25*np.mean(rat_dd_da[-10:])
    # ax.set_xlim(1e-3, 10)
    # ax.set_ylim(0,ymax)
    # ax.set_ylim(1,3)      
    ax.set_ylim(ymin,ymax)
    plt.savefig('FRET-FCS-all.png', dpi=300, transparent=False, bbox_inches='tight')


if 0:
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])

    ncol = len(measurement_group)
    colors = plt.cm.jet(np.linspace(0,1,ncol))
    i = 0
    for m in measurement_group:
        measurement_time = float(re.findall(r'T\d+',m['name'])[0][1:]) #convert to float if want to do maths
        measurement_time = measurement_time/60 #now minutes
        # measurement_time = re.findall(r'T\d+',m['name'])[0][1:] #or keep as string if just for label
        t = m['DxA']['time']
        Gx = m['DxA']['G']
        errx = m['DxA']['err']
        
        t = m['AA']['time']
        Gaa = m['AA']['G']
        erraa = m['AA']['err']
        
        t = m['AA']['time']
        Gdd = m['DD']['G']
        errdd = m['DD']['err']
        
        rat_dd_da = Gdd/Gx
        ax.plot(t,rat_dd_da, linestyle = '', marker = 'o', markersize = 1, label = 't = %.2f min' %(measurement_time), color = colors[i])        
        i+=1        
        
    ax.set_xscale('log')
    # ax.legend()
        
    ax.set_xlabel(r'$\tau$ (ms)', labelpad=10)
    ax.set_ylabel(r'$G(\tau)$', labelpad=10)
    
    ymax = 1.5*np.mean(rat_dd_da[:5])
    ymin = 0.25*np.mean(rat_dd_da[-10:])
    # ax.set_xlim(1e-3, 10)
    # ax.set_ylim(0,4)
    # ax.set_ylim(1,3)      
    ax.set_ylim(ymin,ymax)
    
    plt.savefig('FRET-FCS.png', dpi=300, transparent=False, bbox_inches='tight')