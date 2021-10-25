# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 13:35:51 2021

@author: gwg24
"""
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
sys.path.append("./lib") #point python to location containing the below three modules
import FCS_fitfunc as ff
import SPT_reader as spt
import FCS_helpful as fcs
import lmfit

def Kinetic(t,a,b,c,d,tau):
    
    Xdd = a + b*np.exp(-t/tau) 
    Xda = c + d*np.exp(-t/tau)
    X = Xdd/Xda

    return X

def simplefit(t, X, err):
    # 1-comp diffusion
    
    # A0_guess = np.mean(G[:5])
    # arg = abs(G - A0_guess/2)
    # ind = np.argmin(arg)
    # tau_guess = t[ind]    
    # print('Guessing tau = %f ms' %(tau_guess))
    
    model = lmfit.Model(Kinetic)
    params = model.make_params(a=1, tau=3)
    params['a'].set(value=1)
    params['b'].set(value=1, min = 0)
    params['c'].set(value=1)
    params['d'].set(value=-0.5, max = 0)
    params['tau'].set(min = 0, value = 3)
 
    weights = 1/err

    fitres = model.fit(X, t=t, params=params, method='least_squares',
                       weights=weights)
    print('\nList of fitted parameters for %s: \n' % model.name)
    fitres.params.pretty_print(colwidth=10, columns=['value', 'stderr', 'min', 'max'])
    
    print(fitres.fit_report())
    
    return fitres


t = np.logspace(-5,0,num = 1000)
# t = t *1000
# m = spt.Read_FCS('./Data/SDS_excess_ABO-FRET-FCS.dat')
# m_group = spt.Read_FCS('./Data/Barghorn_4xdil_FRET-FCS.dat')
# m_group = spt.Read_FCS('./Data/SDS_excess_ABO-FRET-FCS.dat')
m_group = spt.Read_FCS('./Data/Barghorn_4xdil_nextday_filt1_FRET-FCS.dat')
# m_group = spt.Read_FCS('./Data/A488_A594_FRET_FCS_grouped.dat')
# m_group = spt.Read_FCS('./Data/Barghorn_4xdil_nextday_4xfilt_FRET-FCS.dat')
# k_21 = 50
# # k_12 = 1/(100E-6)
# k_12 = 20
# E1 = 0.7
# E2 = 0
# Xdd = 1 + ( (k_12*k_21*((E1 - E2)**2) )/( (k_21*(1-E1) + k_12*(1-E2))**2 ) )*np.exp(-(k_12 + k_21)*t)

# Xda = 1 - ( (k_12*k_21*((E1 - E2)**2) )   /( (k_21*(1-E1) + k_12*(1-E2))*(k_21*E1 + k_12*E2)  ))*np.exp(-(k_12 + k_21)*t)

# const = 0.6
# X = Xdd/Xda - const

# plt.figure()
# plt.plot(t*1000,X)

tau = []
for m in m_group[:-1]:
    
    t = m['DD']['time']
    Gdd = m['DD']['G']
    errdd = m['DD']['err']
    
    t = m['DxA']['time']
    Gx = m['DxA']['G']
    errx = m['DxA']['err']
    
    X = Gdd/Gx
    errX = X*np.sqrt( (errdd/Gdd)**2   + (errx/Gx)**2    )
    
    fitres =  simplefit(t, X, errX)
    
    tau.append(fitres.values['tau'])
    
    plt.figure()
    plt.plot(t,X)
    plt.plot(t, fitres.best_fit)
    plt.xscale('log')
    plt.ylim(0.5*np.min(fitres.best_fit),1.5*np.max(fitres.best_fit))
    
    
    
plt.figure()
plt.plot(t, fitres.residual)    
    
tau = np.array(tau)

plt.figure()
plt.plot(tau)
plt.figure()
plt.hist(tau)