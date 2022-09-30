# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 15:29:22 2021

Purpose: To fit autocorrelation curves from SPT (SymPhoTime software)

IN SPT for FRET-FCS
A = Ch2, Al488, 0-50 ns
B = Ch1, Al594, 0-50 ns

In SPT for FCCS
A = Ch2, Al488, 0-50 ns
B = Ch1, Al594, 50-100 ns

Do not export with fits.

Usage: Look for FIXME in code and make changes accordingly to fit
application.

@author: gwg24
"""

import numpy as np
import numpy.ma as ma
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import lmfit
import glob
from uncertainties import ufloat
from uncertainties.umath import *

import sys
sys.path.append("./lib") #point python to location containing the below three modules
import FCS_fitfunc as ff
import platform
if platform.system() == ' Windows':
    import SPT_reader as spt
    print('WINDOWS')  
elif platform.system() == 'Linux':
    import SPT_reader_edit as spt
    print('LINUX')  
elif platform.system() == 'Darwin':
    import SPT_reader_edit as spt
    print('MAC')  
else:
    print('OH NO, YEET, NO OPERATING SYSTEM')  
import FCS_helpful as fcs

# Edit the font, font size, and axes width
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1



'''''''''''''''''''''Functions to facilitate fitting'''''''''''''''''''''''



'''One component model'''

def simplefit(t, G, err, kappa = 5, triplet = False, t_min = 1E-9, t_max = 1e7):
    # 1-comp diffusion, triplet considered only if set to True
    A0_guess = np.mean(G[10:20])
    arg = abs(G - A0_guess/2)
    ind = np.argmin(arg)
    tau_guess = t[ind]    
    print('Guessing tau = %f ms' %(tau_guess))
    
    model = lmfit.Model(ff.diffusion_3d_triplet)
    params = model.make_params(A0=A0_guess, tau_diff=tau_guess)
    #Use triplet model or not.
    if triplet:
        params['T'].set(min=0, max=1, value=0.1)
        params['tau_t'].set(min=1E-9, max=0.05, value=10e-4)
    else:
        params['T'].set(value=0, vary=False)
        params['tau_t'].set(value=1e-4, vary=False)
    params['A0'].set(min=0.00001, value=A0_guess)
    params['tau_diff'].set(min=1e-4, value=tau_guess)
    params['Ginf'].set(value=0, vary = True)
    params['kappa'].set(value=kappa, vary=False)  # 3D model only
    
    weights = 1/err
    #t_max = 1e7
    #FIXME
    fitres = model.fit(G[np.logical_and(t < t_max, t > t_min)], timelag=t[np.logical_and(t < t_max, t > t_min)], params=params, method='least_squares',
                       weights=weights[np.logical_and(t < t_max, t > t_min)])
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
    
    print(fitres.fit_report())
    return fitres

''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''Body'''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''

'''Usage Options. Input required below.'''

## Autocorrelation of donor, receptor, or cross-correlation?
# Possible key values are DD (autocorrelation Donor Donor), AA (auto, accepptor acceptor), DxA (donor acceptor cross correlation)
key = 'DD'

#Fitting Model (Triplet? True/False)
triplet = False

##Min and max for fitting (units of ms).
#Use really low and high values to use all data.
t_min = 1e-4 #ms
t_max = 1e7 #ms

'''Input Paramaters (from calibration-see FCS_calibration.py)'''
set_kappa = 6.7085                   # from calibration, A488
# set_kappa = 6.85                     # from calibration, A594
td_ref = ufloat(0.0381, 0.00018)     # from calibration (ms), A488 
# td_ref = ufloat(0.0659, 0.00077)     # from calibration (ms), A594
D_ref = ufloat(470, 40)              # from literature, for calibration (um^2/s), Rho110, 470 um^2/s, 25 C
# D_ref = ufloat(350, 40)              # from literature, for calibration (um^2/s), a594-SDP, 350 um^2/s, 25 C
temperature_ref = ufloat(25, 0.5)    # temperature at which reference D was taken (celsius)
temperature_lab = ufloat(21.5, 0.5)    # our labs temeprature (celsius)
    
Toffset = 0 #offset the time axis by some amount (make sure consistent with units for measurement_time (i.e., minutes, hours?))

#FIXME, put your own path HERE
path = "/Users/bab226/Documents/yale_research/iapp/fcs/fcs_analyzer_package/fcs_analysis_package/Data/BB_thermo_iao_sds_dilutions_02-11-22.sptw/day1/"
for name in glob.glob(path + "*.dat"):
    name = os.path.basename(name)[:-4]
    measurement_group = spt.Read_FCS(path + name)
    # measurement_group = spt.Read_FCS('')
    # measurement_group = spt.Read_FCS('')
    print("""
    ##############################
    ##### Now Analyzing %s #######
    ##############################
    """ %(name))
    '''Example of batch fitting/plotting'''
    
    #Intialize lists to store results of fitting, 1c = 1 component.
    
    #1c Parameters
    tau_1c = []
    tau_1c_err = []
    tau_trip_1c = []
    tau_trip_1c_err = []
    T_1c = []
    T_1c_err = []
    N_1c = []
    N_1c_err = []
    redchi1c = []
    
    
    t_measurement = [] #measurement time which will be grepped from title of measurement
     
    
    #If plotting time series on same graph, this can be helpful
    #There are n measurements in measurement group, so we divide a color map into n linearly spaced colors
    #i.e., plot goes from hot to cold colors or whatever --can change color map from jet to other
    ncol = len(measurement_group)
    colors = plt.cm.jet(np.linspace(0,1,ncol))
    
    i = 0
    Gsum = 0
    for m in measurement_group[:]: #all of them
    # for m in measurement_group[:-1]: #exclude the last measurement (e.g., incomplete, or too much evaportation or ...)    
    # for m in measurement_group[::2]: #everyother meauserment    
    
        #SPT64 saves file name in measurement group as e.g., 'ThisName_T1800s_1.ptu'
        #where T1800s means this measurement was taken at 1800 seconds after start.
        #Use regular expression to find pattern Tsomeinteger and then convert it to a number
        if len(measurement_group) == 1:
            measurement_time = 0
        else:
            measurement_time = float(re.findall(r'T\d+',m['name'])[0][1:])
            measurement_time = Toffset + (measurement_time/60) #now minutes 
        mylabel = 't = %.2f min' %(measurement_time)
        
        #Use key to decide which correlation curve to look at    
        t = m[key]['time']
        G = m[key]['G']
        err = m[key]['err']
        
        Gsum = Gsum + G    
        
        print()
        print()
        print('############')
        print('Results of measurement at time t = %.2f min' %(measurement_time))
        fitres1c = simplefit(t, G, err, kappa = set_kappa, triplet = triplet, t_min = t_min, t_max = t_max)

        redchi1c.append(fitres1c.redchi)
        # redchi2c.append(fitres2c.redchi)
        
        #Store results
        t_measurement.append(measurement_time)
        
        tau_1c.append(fitres1c.values['tau_diff'])
        tau_1c_err.append(fitres1c.params['tau_diff'].stderr)
        tau_trip_1c.append(fitres1c.values['tau_t'])
        tau_trip_1c_err.append(fitres1c.params['tau_t'].stderr)
        T_1c.append(fitres1c.values['T'])
        T_1c_err.append(fitres1c.params['T'].stderr)
        A01c = ufloat(fitres1c.values['A0'], fitres1c.params['A0'].stderr) 
        x = 1/A01c #N = 1/A
        N_1c.append(x.nominal_value)
        N_1c_err.append(x.std_dev)
     
        width = 3.42
        fig, ax = plt.subplots(2, 1, figsize=(width,width/1.2), sharex=True,
                           gridspec_kw={'height_ratios': [3, 1.2]})
        plt.subplots_adjust(hspace=0.3)
        ax[0].errorbar(t,G, yerr = err, linewidth =1, label = mylabel, linestyle = '', marker = 'o', markersize = 2, color = colors[i])
        
        
        ax[0].plot(t[np.logical_and(t < t_max, t > t_min)], fitres1c.best_fit, color = 'm', label = '1-comp', linestyle = '--', zorder = 10) #zorder forces best fit to plot on top errorbar
    
        # resid1 = (G - fitres1c.best_fit)/err
        # resid2 = (G - fitres2c.best_fit)/err
        
        ax[0].set_xscale('log')
        ax[1].plot(t[np.logical_and(t < t_max, t > t_min)], fitres1c.residual, 'm')
       
        # ax[1].plot(t, resid1, 'm')
        
        ax[1].set_xscale('log')
        ax[0].legend()
       
        std_wres1 = np.std(fitres1c.residual)
      
        # ax[1].set_ylim(-ym, ym)
        max_std = np.std(fitres1c.residual)
        ax[1].set_ylim(- 3*max_std, + 3*max_std)
        ax[1].set_ylim(-5,5)
        ax[1].set_xscale('log')
        
        for a in ax:
            a.grid(True); a.grid(True, which='minor', lw=0.3)
        ax[1].set_xlabel(r'$\tau$ (ms)', labelpad=5)
        ax[0].set_ylabel(r'$G(\tau)$', labelpad=5)
        # ax[0].set_ylabel('G(τ)')
        ax[1].set_ylabel('wres', labelpad = 5)
        # ax[0].set_title('Pseudo Autocorrelation')
        # ax[1].set_xlabel('Time Lag, τ (s)');
        
        savefig_string = './Figures/' + m['name'][:-4] + '_' + key + '.png' #saves in Figures folder, as measurement name, appended with which correlation
        plt.savefig(savefig_string, dpi=300, transparent=False, bbox_inches='tight')
        i += 1
    
    Gsum = Gsum/len(measurement_group)
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.plot(t,Gsum, linewidth =1, linestyle = '', marker = 'o', markersize = 2, color = 'k')
    ax.set_xscale('log')
    plt.title('Average correlation curve over all measurements')

    #turn results into np arrays so we can do maths with them
    t_measurement = np.array(t_measurement)
    
    #Veff calculation
    w = (4*D_ref*(td_ref/1000))**(0.5) #um   
    w = w*(1e-6) #now in meters
    #print(w)
    kappa = set_kappa
    Veff = ((np.pi)**(3/2)) * kappa * (w**3) #in m^3
    Veff = Veff *1000 #now in L 
    
    N = ufloat(np.mean(N_1c), np.std(N_1c)) #particles
    N = N/(6.022E23) #particles/(particles/mol)  = mol

    C = N/Veff #mol / L = M
    
    #See NEW function in FCS_helpful.py
    #Veff = fcs.get_veff(td_ref, D_ref, set_kappa)  #in L #Greg, check maths, use real
    
    print('Veff in fL')
    print(Veff*1E15)
    
    #Calculate concentration from N_2c
    Na = 6.022*10**23
    
    #Chi-squared from least-squares fitting 
    redchi1c = np.array(redchi1c)
    
    tau_1c = np.array(tau_1c)
    tau_1c_err = np.array(tau_1c_err)
    N_1c = np.array(N_1c)
    N_1c_err = np.array(N_1c_err)
    N_1c = np.array(N_1c)
    N_1c_err = np.array(N_1c_err)
    
    ##NEED ATTENTION for same reasons as above
    #Calculate concentration from N_1c
    Na = 6.022*10**23
    N_mol = N_1c/Na
    N_mol_err = N_1c_err/Na
    C_1c = N_mol/Veff.nominal_value*1E9 #nM
    C_1c_err = N_mol_err/Veff.nominal_value*1E9 #nM
    
    D1c, eD1c = fcs.td2D(tau_1c, tau_1c_err, temperature_lab = temperature_lab, td_ref = td_ref, D_ref = D_ref, temperature_ref = temperature_ref )
    Rh1c, err_Rh1c = fcs.D2Rh(D1c, eD1c, temperature_lab = temperature_lab)
    
    print()
    print()
    print('############')
    print('Quick summary')
    print('Analyzing ' + key + ' curves')
    print('Mean 1-comp diffusion time = %.4f +/- %.5f ms' %(np.mean(tau_1c), np.std(tau_1c, ddof =1)))
    print('Mean 1-comp diffusion coeff = %.4f +/- %.5f um^2 s^-1' %(np.mean(D1c), np.std(D1c, ddof =1)))
    print('Mean 1-comp Rh = %.4f +/- %.5f nm' %(np.mean(Rh1c), np.std(Rh1c, ddof =1)))
    print('Mean 1-comp N = %.4f +/- %.5f nm' %(np.mean(N_1c), np.std(N_1c, ddof =1)))
    print('Mean 1-comp redchi = %.4f +/- %.5f' %(np.mean(redchi1c), np.std(redchi1c, ddof =1)))
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
    ax.errorbar(t_measurement, tau_1c*1000, yerr = tau_1c_err*1000, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
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
    ax.errorbar(t_measurement, T_1c, yerr = T_1c_err, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'Fraction trip', labelpad=10)
    plt.savefig('./Figures/frac_trip_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    width = 3.42
    fig = plt.figure(figsize=(width,width/1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.errorbar(t_measurement, tau_trip_1c, yerr = tau_trip_1c_err, linewidth =1, label = '', linestyle = '', marker = 'o', markersize = 4, color = color_dict[key])
    ax.set_xlabel(r'time (min)', labelpad=10)
    ax.set_ylabel(r'$t_{trip}$ ($\mathrm{\mu s}$)', labelpad=10)
    plt.savefig('./Figures/tau_trip_' + key + '.png', dpi=300, transparent=False, bbox_inches='tight')
    
    
    
    # #Diffusion Time 1
    # with open('./Results/'+ name + '_tD_fast_' + key + '.dat', "w" ) as f:
    #     f.write('t \t T1 \t err \n')
    #     for t,R,err in zip(t_measurement, tau_fast2c, tau_fast_err2c):
    #         f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))
     
        
    # #Diffusion Time 2
    # with open('./Results/'+ name + '_tD_slow_' + key + '.dat', "w" ) as f:
    #     f.write('t \t T2 \t err \n')
    #     for t,R,err in zip(t_measurement, tau_slow2c, tau_slow_err2c):
    #         f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))
        
    # # #Diffusion Coefficient 1
    # with open('./Results/'+ name + '_D2c_fast_' + key + '.dat', "w" ) as f:
    #     f.write('t \t D1 \t err \n')
    #     for t,R,err in zip(t_measurement, D2c_fast, eD1c):
    #         f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))
            
    # # #Diffusion Coefficient 2
    # with open('./Results/'+ name + '_D2c_slow_' + key + '.dat', "w" ) as f:
    #     f.write('t \t D2 \t err \n')
    #     for t,R,err in zip(t_measurement, D2c, eD2c):
    #         f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))
            
    # # #Hydrodynamic radius
    # with open('./Results/'+ name + '_Rh2c_' + key + '.dat', "w" ) as f:
    #     f.write('t \t R \t err \n')
    #     for t,R,err in zip(t_measurement, Rh2c,err_Rh2c):
    #         f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))
              
    # ############
    # ##Averages##
    # ############
    
    # # Sample Name: %s
    # # Mean 2-comp fast diffusion time = %.4f +/- %.5f ms
    # # Mean 2-comp slow diffusion time = %.4f +/- %.5f ms
    # # Mean 2-comp fast diffusion coeff = %.4f +/- %.5f um^2 s^-1
    # # Mean 2-comp slow diffusion coeff = %.4f +/- %.5f um^2 s^-1
    # # Mean 2-comp Rh slow = %.4f +/- %.5f nm
    # # Mean 2-comp f(slow) = %.3f +/- %.4f
    # # Mean 2-comp (chi-squared) redchi = %.4f +/- %.5f 
    # # """ %(name, np.mean(tau_fast2c), np.std(tau_fast2c, ddof=1), np.mean(tau_slow2c),
    #     np.std(tau_slow2c, ddof =1), np.mean(D1c), np.std(D1c, ddof =1), 
    #     np.mean(D2c), np.std(D2c, ddof =1), np.mean(Rh2c), np.std(Rh2c, ddof =1), 
    #     np.mean(p1), np.std(p1, ddof=1), np.mean(redchi2c), np.std(redchi2c, ddof =1)),
    #     file=f)
    
    
    #Diffusion Time
    with open('./Results/' + name + '_Dt_' + key + '.dat', "w" ) as f:
        f.write('t \t T \t err \n')
        for t,R,err in zip(t_measurement, tau_1c, tau_1c_err):
            f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))
       
      
    #Diffusion Coefficient
    with open('./Results/' + name + '_D1c_' + key + '.dat', "w" ) as f:
        f.write('t \t D \t err \n')
        for t,R,err in zip(t_measurement, D1c, eD1c):
            f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))  
       
       
    #Hydrodynamic Radius
    with open('./Results/' + name + '_Rh1c_' + key + '.dat', "w" ) as f:
        f.write('t \t R \t err \n')
        for t,R,err in zip(t_measurement, Rh1c, err_Rh1c):
            f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))    
       
    #Number of Molecules
    with open('./Results/' + name + '_N1c_' + key + '.dat', "w" ) as f:
        f.write('t \t N \t err \n')
        for t,R,err in zip(t_measurement, N_1c, N_1c_err):
            f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))    
            
    #Concentration
    with open('./Results/' + name + '_C1c_' + key + '.dat', "w" ) as f:
        f.write('t \t C \t err \n')
        for t,R,err in zip(t_measurement, C_1c, C_1c_err):
            f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))  
            
    #Chi-squared 1c
    with open('./Results/' + name + '_redchi1c_' + key + '.dat', "w" ) as f:
        f.write('t \t redchi \n')
        for t,R in zip(t_measurement, redchi1c):
            f.write('%.3f \t %.3f \n' %(t,R)) 
            
    with open('./Results/' + name + '_summary1c_' + key + '.dat', "w" ) as f:
        
        print("""
    ############
    ##Averages##
    ############
    
    Sample Name: %s
    Mean 1-comp diffusion time = %.4f +/- %.5f ms
    Mean 1-comp diffusion coeff = %.4f +/- %.5f um^2 s^-1
    Mean 1-comp Rh = %.4f +/- %.5f nm
    Mean 1-comp redchi = %.4f +/- %.5f #Chi squared statistic
    Mean 1-comp number of molecules = %.3f +/- %.4f
    Mean 1-comp sample conc (nM) = %.3f +/- %.3f
    """ %(name, np.mean(tau_1c), np.std(tau_1c, ddof =1), np.mean(D1c), 
            np.std(D1c, ddof =1), np.mean(Rh1c), np.std(Rh1c, ddof =1), 
            np.mean(redchi1c), np.std(redchi1c, ddof =1), np.mean(N_1c), 
            np.std(N_1c, ddof = 1), np.mean(C_1c), np.std(C_1c, ddof = 1)), 
            file=f)
    
        
    ''' For MaxEnt '''
    with open('./Results/' + name + '_MaxEnt_' + key + '.dat', "w" ) as f:
        f.write('w (nm) \t trip (us) \t T \t kappa \t D_1c \t err \n')
        for i in range(0,len(tau_trip_1c)):
            f.write('%.2f \t %.2f \t %.3f \t %.2f \t %.2f \t %.2f \n' %(w.nominal_value*10**9,tau_trip_1c[i]*1000, T_1c[i], kappa, D1c[i], eD1c[i]))
        f.write('%.2f \t %.2f \t %.3f \t %.2f \t %.2f \t %.2f \n' %(w.nominal_value*10**9,np.mean(tau_trip_1c)*1000, np.mean(T_1c), kappa, np.mean(D1c), np.mean(eD1c)))
        