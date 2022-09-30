# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 15:38:11 2021

@author: gwg24
"""

import numpy as np
from math import e
# import matplotlib.pyplot as plt
# from uncertainties import ufloat
# from uncertainties.umath import *

# import scipy as sp

''' Basic Diffusion models'''
def diffusion_3d(timelag, tau_diff, A0, Ginf, kappa):
    #1component diffusion 3-D
       return ( Ginf + A0 * 1/(1 + timelag/tau_diff) *
               1/np.sqrt(1 + timelag/(tau_diff*(kappa**2))))
 
def twocomp_diff3d(timelag, p1, tau_diff1, tau_diff2, A0, Ginf, kappa):    
    #Two component diffusion 3-D
    return (A0*(diffusion_3d(timelag, tau_diff1, p1, 0, kappa) + diffusion_3d(timelag, tau_diff2, 1 - p1, 0, kappa)) + Ginf)    

''' Triplet and Diffusion model '''

def triplet(timelag, tau_t, T):
    return (1- T + T*np.exp(-timelag/tau_t))*(1/(1-T))

def diffusion_3d_triplet(timelag, tau_diff, tau_t, T, A0, Ginf, kappa):
    return (diffusion_3d(timelag, tau_diff, A0, Ginf, kappa) * triplet(timelag, tau_t, T))

def twocomp_diffusion_3d_triplet(timelag, p1, tau_diff1, tau_diff2, tau_t, T, A0, Ginf, kappa):
    return (twocomp_diff3d(timelag, p1, tau_diff1, tau_diff2, A0, Ginf, kappa) * triplet(timelag, tau_t, T))


''' Kinetic contribution due to FRET'''
def E_dd(timelag, tau_k, a, b ):
    
    return (1 + a*np.exp(-timelag/tau_k) + b)

def E_aa(timelag, tau_k, c, d ):
    
    return (1 + c*np.exp(-timelag/tau_k) + d)

def E_cross(timelag, tau_k, a, c, b, d):
    
    return (1 - np.sqrt(a*c)*np.exp(-timelag/tau_k) - np.sqrt(b*d))

''' Total correlation models 1-comp FRET+Diffusion'''
def Gdd_model(timelag,  tau_diff, A0, Ginf, kappa, tau_k, a, b):

    return diffusion_3d(timelag, tau_diff, A0, Ginf, kappa)*E_dd(timelag, tau_k, a, b )

def Gaa_model(timelag,  tau_diff, A0, Ginf, kappa, tau_k, c, d):

    return diffusion_3d(timelag, tau_diff, A0, Ginf, kappa)*E_aa(timelag, tau_k, c, d )    

def Gx_model(timelag, tau_diff, A0, Ginf, kappa, tau_k, a, c, b, d):
    return diffusion_3d(timelag, tau_diff, A0, Ginf, kappa)*E_cross(timelag, tau_k, a, c, b, d)

''' Total correlation models 2-comp FRET+Diffusion'''

def Gdd_model2(timelag,  p1, tau_diff1, tau_diff2, A0, Ginf, kappa, tau_k, a, b):

    return twocomp_diff3d(timelag, p1, tau_diff1, tau_diff2, A0, Ginf, kappa)*E_dd(timelag, tau_k, a, b )

def Gaa_model2(timelag,  p1, tau_diff1, tau_diff2, A0, Ginf, kappa, tau_k, c, d):

    return twocomp_diff3d(timelag, p1, tau_diff1, tau_diff2, A0, Ginf, kappa)*E_aa(timelag, tau_k, c, d )

def Gx_model2(timelag,  p1, tau_diff1, tau_diff2, A0, Ginf, kappa, tau_k, a,c,b,d):

    return twocomp_diff3d(timelag, p1, tau_diff1, tau_diff2, A0, Ginf, kappa)*E_cross(timelag, tau_k, a,c, b,d )

''' Global fitting functions FRET-FCS'''

def global_E_model(timelag, tau_diff, A0, Ginf, kappa, tau_k, a, c, b, d):
    #1-comp diffusion + FRET-FCS
    Gdd = Gdd_model(timelag,  tau_diff, A0, Ginf, kappa, tau_k, a, b)
    Gaa = Gaa_model(timelag,  tau_diff, A0, Ginf, kappa, tau_k, c, d)
    Gx = Gx_model(timelag, tau_diff, A0, Ginf, kappa, tau_k, a, c, b, d)
    
    #do global fit by concatenating residuals, i.e., data and model are concatenation of the three Gs w/ shared parameter
    G_global = np.concatenate((Gdd,Gaa,Gx), axis = None)
    return G_global

def global_E_model2(timelag,  p1, tau_diff1, tau_diff2, A0, Ginf, kappa, tau_k, a, c, b, d):
    #two-comp diffusion + FRET-FCS
    Gdd = Gdd_model2(timelag,  p1, tau_diff1, tau_diff2, A0, Ginf, kappa, tau_k, a, b)
    Gaa = Gaa_model2(timelag,  p1, tau_diff1, tau_diff2, A0, Ginf, kappa, tau_k, c, d)
    Gx = Gx_model2(timelag,  p1, tau_diff1, tau_diff2, A0, Ginf, kappa, tau_k, a,c,b,d)
    
    #do global fit by concatenating residuals, i.e., data and model are concatenation of the three Gs w/ shared parameter
    G_global = np.concatenate((Gdd,Gaa,Gx), axis = None)
    return G_global

def global_E_data(Gdd,Gaa,Gx, errdd, erraa, errx):    
    #prepare vector of data in same order as model
    G_global = np.concatenate((Gdd,Gaa,Gx), axis = None)
    
    err_global = np.concatenate((errdd,erraa,errx), axis = None)
    return G_global, err_global

''' Guassian fit for maximum entropy analysis '''
def gaussian(log_D, a, mu, sigma):
    #Single gaussian form
    return (a*np.e**((-(log_D-mu)**2)/2*sigma**2))
