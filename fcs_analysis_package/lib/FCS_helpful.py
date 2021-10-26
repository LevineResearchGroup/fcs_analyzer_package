# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 15:05:31 2021

@author: gwg24 (Greg Gomes)

Updated with more functionalities by Bryan Bogin.


Important note: viscosity here is temperature corrected water. High concentrations NaCl, GdmCl, glycerol, crowding agent etc., would violate that.

"""


import numpy as np

from uncertainties import ufloat
from uncertainties.umath import *

def td2D(td_array, etd_array,temperature_lab, td_ref, D_ref, temperature_ref = 25 ):
    
    #First need to correct the reference D to lab temperature
    T_lab = temperature_lab + 273.15
    T_ref = temperature_ref + 273.15
    A=2.414E-5 # Pa*S
    B=247.8 #
    C=140;
    eta_lab=A*10**(B/(T_lab-C)) # Pa*s
    eta_ref=A*10**(B/(T_ref-C)) # Pa*s
    # k = 1.3807E-23 #Boltzman constant J/K
    D_ref_lab = D_ref * (T_lab/ eta_lab) * (eta_ref / T_ref) #diffusion coeffecient of reference compoound at actual lab temperature
    print(D_ref_lab)
    D_array = []
    eD_array = []
    for td, etd in zip(td_array, etd_array):
        td_sample = ufloat(td,etd)
        D_sample =  D_ref_lab *td_ref/td_sample
        D_array.append(D_sample.nominal_value)
        eD_array.append(D_sample.std_dev)
    
    return np.array(D_array), np.array(eD_array)
    

def D2Rh(D_array,eD_array, temperature_lab = 22):
    #D should be in micrometers^2 / s
    #temperature in degrees C
    #Assumes eta = viscocity of water, corrected for temperature
    
    #Output Rh in nm
    
    Rh_array = []
    Rh_err_array = []
    for D,eD in zip(D_array, eD_array):
        
        D = ufloat(D,eD)
        D = D*(1E-6)**2 #now m^2/s 
                
        T = temperature_lab +273.15
        A=2.414E-5 # Pa*S
        B=247.8 #
        C=140;
        eta=A*10**(B/(T-C)) # Pa*s
        k = 1.3807E-23 #Boltzman constant J/K
    
        Rh = k*T / (6*np.pi*eta*D) #meters
        
        Rh_array.append(Rh.nominal_value*1E9)
        Rh_err_array.append(Rh.std_dev*1E9)
    
    return np.array(Rh_array), np.array(Rh_err_array)

def set_wavelength(laser):
    '''Options are "blue", "green", or "both"'''
    laser = "both"
    if laser == "blue":
        wavelength = 485
    if laser == "green":
        wavelength = 599
    if laser == "both":
        wavelength = np.sqrt(2)*((485*599)/(np.sqrt(485**2+599**2)))
    
    return wavelength

def get_veff(td_ref, D_ref, set_kappa):
    '''Returns Veff in L "'''
    w = (4*D_ref*(td_ref/1000))**(0.5) #um   
    w = w*(1e-6) #now in meters
    #print(w)
    kappa = set_kappa
    Veff = ((np.pi)**(3/2)) * kappa * (w**3) #in m^3
    Veff = Veff *1000 #now in L 
    
    return Veff



