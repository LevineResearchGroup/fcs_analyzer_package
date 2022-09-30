# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 14:35:29 2021

@author: gwg24
"""

import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# from uncertainties import ufloat
# from uncertainties.umath import *
import pandas as pd
import seaborn as sns
from scipy.ndimage.filters import uniform_filter1d

mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
plt.rcParams["errorbar.capsize"] = 4
lw = 2

colors = sns.color_palette("Set2", 5 )
# colors = sns.color_palette("rocket", 3 )
color1 = colors[0]
color2 = colors[1]
color3 = colors[2]

color1 = (27/255,158/255,119/255)
color2 = (217/255, 95/255, 2/255)
color3 = (117/255, 112/255, 179/255)


colorlist = sns.color_palette("Blues_r", n_colors = 5)
# colorlist = colors

# color2 = 'palevioletred'
# color1 = 'lightskyblue'


def ReadData(filename):
    
    '''
    
    '''
    
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


    # print(R)
    return np.array(t), np.array(R), np.array(err)
    # return np.array(t)[:15], np.array(R)[:15], np.array(err)[:15]

def ReadData2(filename):
    
    '''
    WARNING: I ADDED A LIMITER (ONLY RETURN FIRST FIVE) BECAUSE I WANTED EACH TO HAVE ONLY 5 DATA POINTS!
    '''
    
    t = []
    redchi = []
    # err = []
    hl=1
    with open(filename, newline = '') as f:                                                                                          
        reader = csv.reader(f, delimiter='\t')
        for j in range(hl):
            #just skip headerlines, already stored
            next(reader)
        for row in reader:
            #Read remaining lines (row at a time)
            
            t.append(float(row[0]))
            redchi.append(float(row[1]))


    # print(R)
    return np.array(t), np.array(redchi)


plt.rc('font', family='serif',serif='Times')
# plt.rc('font', family='serif')
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes',labelsize=20)


key = 'AA'
path = './Results/iapp_pbs_kinetics_fccs/'
name = 'iapp_pbs_kinetics'

''' data set 1'''
t1, D1, errD1 = ReadData(path + name + '_D1c_' + key + '.dat')
t1, R1, errR1 = ReadData(path + name + '_Rh1c_' + key + '.dat')
t1, N1, errN1 = ReadData(path + name + '_N1c_' + key + '.dat')
t1, redchi1_1 = ReadData2(path + name + '_redchi1c_'+ key +'.dat')
#t2, D2, errD2 = ReadData(path + name + '_D2c_slow_' + key + '.dat')
t2, R2, errR2 = ReadData(path + name + '_Rh2c_' + key + '.dat')
t2, N2, errN2 = ReadData(path + name + '_N2c_' + key + '.dat')
#t2, P1_2, errP1_2 = ReadData(path + name + '_p1_2c_' + key + '.dat')    #Co-diffusion fraction

# t1, p1, p1_err = ReadData(path + name + '_p1_' + key + '.dat')

'''data set 2, to be joined with 1'''
# t1_2, R1_2, err1 = ReadData('../Results/iapp_pbs_kinetics_fccs/iapp_pbs_kinetics_Rh1c_'+ key+'.dat')
# t1_2, redchi1_2 = ReadData2('./Results/ApoE2_1uL_1mL_PBS_group2_redchi1c_'+ key+'.dat')
# t1_2, N1_2, errN1 = ReadData('./Results/ApoE2_1uL_1mL_PBS_group2_N_1c_' + key + '.dat')

''' manually putting in keys for cross correlation'''
# t1, N1_x, errN1 = ReadData('./Results/iapp_pbs_kinetics_fccs/iapp_pbs_kinetics_N_1c_' + 'DxA' + '.dat')
#t1, N1_D, errN1 = ReadData(path + name + '_N1c_' + 'DD' + '.dat')
#t1, N1_A, errN1 = ReadData('./Results/iapp_pbs_kinetics_fccs/iapp_pbs_kinetics_N_1c_' + 'AA' + '.dat')


# print(len(t1))
# exit()

# plt.figure()
# fig, ax = plt.subplots()
# fig.subplots_adjust(left=0, bottom=.15, right=.99, top=.97,hspace = 0, wspace = 0)

'''defintion of fraction codiffusing'''
# fD = (1/N1_x)*N1_D
# fA = (1/N1_x)*N1_A

# plt.plot(t1,fD, 'bo', alpha = 0.2, label = 'fraction D codiff')
# plt.plot(t1,uniform_filter1d(fD, size  = 10), color = 'b')
# plt.plot(t1,fA, 'go',alpha = 0.2, label = 'fraction A codiff')
# plt.plot(t1,uniform_filter1d(fA, size  = 10), color = 'g')
# plt.xlabel('time (hrs)')
# plt.ylabel(r'$f_{x}$ ')
# plt.legend()
# width = 8
# height = width /1.2

# fig.set_size_inches(width, height)
# fig.savefig('ApoE-Rh-meas-'+ 'fcodiff' +'.png', dpi = 300, bbox_inches='tight')


#I had two experiments that I needed to join together
# t1_2 = t1_2 + max(t1)
# t = np.concatenate([t1,t1_2])
# R = np.concatenate([R1,R1_2])
# N = np.concatenate([N1,N1_2])



#Calibration/Control
# t1_R110, D1_R110, err_D1_R110 = ReadData('./Results/rho110_control/rho110_D1c_DD.dat')
# t1_R110, R1_R110, err_R1_R110 = ReadData('./Results/rho110_control/rho110_Rh1c_DD.dat')
# t1_R110, N1_R110, err_N1_R110 = ReadData('./Results/rho110_control/rho110_N1c_DD.dat')
# t1_R110, redchi1_R110 = ReadData2('./Results/rho110_control/rho110_redchi1c_DD.dat')

# N1_R110 = N1_R110/(6.022E23) #particles/(particles/mol)  = mol
# N1_R110 = N1_R110/Veff #mol / L = M
# N1_R110 = N1_R110*1e9 #now in nM

# plt.figure()
# y1 = uniform_filter1d(D1_R110, size  = 5, mode='reflect')
# plt.plot(t1_R110,D1_R110,'o', alpha = 0.1, color = 'k')
# plt.plot(t1_R110, y1, '-', alpha = 0.5, color='b' )

# plt.figure()
# y2 = uniform_filter1d(N1_R110, size  = 5)
# plt.plot(t1_R110,N1_R110,'o', alpha = 0.1, color = 'r')
# plt.plot(t1_R110,y2, color = 'r')

# plt.figure()
# y = uniform_filter1d(R1_R110, size = 5)
# plt.plot(t1_R110,R1_R110,'o', alpha = 0.1, color = 'blue')
# plt.plot(t1_R110,y, color = 'b')


# t = np.concatenate([t1,t1_2])
# R = np.concatenate([R1,R1_2])
# redchi = np.concatenate([redchi1_1,redchi1_2])






######## ONE component fitting #########
#Truncate data set here (optional)
t1 = t1
t1 = t1/60 #now in hours
R1 = R1
N1 = N1
D1 = D1
redchi1c = redchi1_1

Veff = 0.62E-15 #V in L
N1 = N1/(6.022E23) #particles/(particles/mol)  = mol
N1 = N1/Veff #mol / L = M
N1 = N1*1e9 #now in nM

#Running average
k = 10
# filt = redchi1c < k
# D1 = D1[filt]
# R1 = R1[filt]
# t1 = t1[filt]
# N1 = N1[filt]
redchi1c_filt = redchi1c

#Plot Diffusion data
# plt.figure()

# fig, ax = plt.subplots()
# fig.subplots_adjust(left=0, bottom=.15, right=.99, top=.97,hspace = 0, wspace = 0)

# y_avg = uniform_filter1d(D1, size  = 50)
# plt.plot(t1,D1,'o', alpha = 0.1, color = 'b', label = 'D1c')
# plt.plot(t1, y_avg, color = 'k', label = 'running average k = %d' %(k))
# #plt.plot(t, y, color = 'b', label = r'running average k = %d, $\chi^2 < %.1f$ ' %(k_filt,thr))
# #plt.ylim(200,350)
# plt.xlabel('time (hrs)')
# plt.ylabel(r'$D$ ($\mathrm{\mu m^2 s^{-1}}$')
# plt.legend()
# print(np.mean(D1))

# width = 8
# height = width /1.2

# fig.set_size_inches(width, height)
#fig.savefig(path+name+'iapp_bps_kinetics_D1c_'+ key+'.png', dpi = 300, bbox_inches='tight')

#Plot Rh Data
plt.figure()
y_avg = uniform_filter1d(R1, size  = k)
plt.plot(t1,R1,'o', alpha = 0.1, color = 'r', label = 'Rh (1-component fitting)')
plt.plot(t1, y_avg, color = 'k', label = 'running average k = %d' %(k))
plt.ylim(0.5,150)
print(np.mean(R1))

filt = redchi1c < 100

R1_filt = R1[filt]
t1_filt = t1[filt]

y_avg_filt = uniform_filter1d(R1_filt, size  = k)
plt.plot(t1_filt,R1_filt,'o', alpha = 0.1, color = 'b', label = 'Filtered Rh (1-component fitting)')
plt.plot(t1_filt, y_avg_filt, color = 'g', label = 'running average k = %d' %(k))

plt.xlabel('time (hrs)')
plt.ylabel(r'$R_{h}$ (nm)')
plt.legend()
print(np.mean(R1_filt))
plt.savefig(path+name+'_Rh_'+ key+'.png', dpi = 300, bbox_inches='tight')

# #Plot # molecules
# plt.figure()
# y_avg = uniform_filter1d(N1, size  = 10)
# plt.plot(t1,N1,'o', alpha = 0.1, color = 'r', label = 'R1c')
# plt.plot(t1, y_avg, color = 'k', label = 'running average k = %d' %(k))
# #plt.ylim(0.5,1)
# plt.xlabel('time (hrs)')
# plt.ylabel('C (nM)')
# plt.legend()
# print(np.mean(N1))

# #Plot frac co-diffusing
# plt.figure()
# y_avg = uniform_filter1d(P1_2, size  = 10)
# plt.plot(t1,P1_2,'o', alpha = 0.1, color = 'r', label = 'frac co-diffusing')
# plt.plot(t1, y_avg, color = 'k', label = 'running average k = %d' %(k))
# #plt.ylim(0.5,1)
# plt.xlabel('time (hrs)')
# plt.ylabel('Co-diffusing fraction')
# plt.legend()
# print(np.mean(N1))

######## Two component fitting #########
#Filter for worst fits

# t2 = t2
# t2 = t2/60 #now in hours
# #D2 = D2
# R2 = R2
# N2 = N2
# #P2 = P2
# redchi1c = redchi1_1

# filt = R2 < 500
# #filt = redchi1c > k
# # D2 = D2[filt]
# R2 = R2[filt]
# t2 = t2[filt]
# N2 = N2[filt]
# P2 = P2[filt]
# redchi2c = redchi1c[filt]

#Plot Diffusion data
# plt.figure()

# fig, ax = plt.subplots()
# fig.subplots_adjust(left=0, bottom=.15, right=.99, top=.97,hspace = 0, wspace = 0)

# y_avg = uniform_filter1d(D2, size  = 30)
# plt.plot(t2,D2,'o', alpha = 0.1, color = 'b', label = 'D2c')
# plt.plot(t2, y_avg, color = 'k', label = 'running average k = %d' %(k))
# #plt.plot(t, y, color = 'b', label = r'running average k = %d, $\chi^2 < %.1f$ ' %(k_filt,thr))
# #plt.ylim(200,350)
# plt.xlabel('time (hrs)')
# plt.ylabel(r'$D$ ($\mathrm{\mu m^2 s^{-1}}$')
# plt.legend()
# print(np.mean(D2))

# width = 8
# height = width /1.2

# fig.set_size_inches(width, height)
#fig.savefig(path+name+'iapp_bps_kinetics_D1c_'+ key+'.png', dpi = 300, bbox_inches='tight')

#Plot Rh Data
# plt.figure()
# y_avg = uniform_filter1d(R2, size  = 10)
# plt.plot(t2,R2,'o', alpha = 0.1, color = 'r', label = 'R2c')
# plt.plot(t2, y_avg, color = 'k', label = 'running average k = %d' %(k))
# plt.ylim(0.5,200)
# plt.xlabel('time (hrs)')
# plt.ylabel(r'$R_{h}$ (nm)')
# plt.legend()
# print(np.mean(R2))

# #Plot # molecules
# plt.figure()
# y_avg = uniform_filter1d(N2, size  = 10)
# plt.plot(t2,N2,'o', alpha = 0.1, color = 'r', label = 'N2c')
# plt.plot(t2, y_avg, color = 'k', label = 'running average k = %d' %(k))
# #plt.ylim(0.5,1)
# plt.xlabel('time (hrs)')
# plt.ylabel('C (nM)')
# plt.legend()
# print(np.mean(N2))

#Plot frac slow molecules
# plt.figure()
# y_avg = uniform_filter1d(P2, size  = 30)
# plt.plot(t2,P2,'o', alpha = 0.1, color = 'r', label = 'Frac slow')
# plt.plot(t2, y_avg, color = 'k', label = 'running average k = %d' %(k))
# #plt.ylim(0.5,1)
# plt.xlabel('time (hrs)')
# plt.ylabel('C (nM)')
# plt.legend()
# print(np.mean(P2))

# #Plotting chi-squared with Running average:
# k = 10
# #k_filt = 5

# plt.figure()

# fig, ax = plt.subplots()
# fig.subplots_adjust(left=0, bottom=.15, right=.99, top=.97,hspace = 0, wspace = 0)

# y = uniform_filter1d(redchi_filt, size  = 10)
# #y_filt = uniform_filter1d(R, size  = 5)
# plt.plot(t_filt, redchi_filt,'o', alpha = 0.1, color = 'b', label = 'P(slow)')
# #plt.plot(t,R,'o', alpha = 0.2, color = 'b', label = 'R')
# plt.plot(t_filt, y, color = 'k', label = 'running average k = %d' %(k))
# #plt.plot(t, y, color = 'b', label = r'running average k = %d, $\chi^2 < %.1f$ ' %(k_filt,thr))
# plt.ylim(0,20)
# plt.xlabel('time (hrs)')
# plt.ylabel('redchi_1c')
# plt.legend()
# print(np.mean(redchi_filt))

# width = 8
# height = width /1.2

# fig.set_size_inches(width, height)
#fig.savefig(path+name+'iapp_bps_kinetics_p1_'+ key+'.png', dpi = 300, bbox_inches='tight')

# #Plotting concentration with Running average:
# k = 10
# #k_filt = 5

# plt.figure()

# fig, ax = plt.subplots()
# fig.subplots_adjust(left=0, bottom=.15, right=.99, top=.97,hspace = 0, wspace = 0)

# y = uniform_filter1d(N_filt, size  = 10)
# #y_filt = uniform_filter1d(R, size  = 5)
# plt.plot(t_filt, N_filt,'o', alpha = 0.1, color = 'b', label = 'P(slow)')
# #plt.plot(t,R,'o', alpha = 0.2, color = 'b', label = 'R')
# plt.plot(t_filt, y, color = 'k', label = 'running average k = %d' %(k))
# #plt.plot(t, y, color = 'b', label = r'running average k = %d, $\chi^2 < %.1f$ ' %(k_filt,thr))
# plt.ylim(4.5,6)
# plt.xlabel('time (hrs)')
# plt.ylabel('N (nM)')
# plt.legend()
# print(np.mean(N_filt))

# width = 8
# height = width /1.2

# fig.set_size_inches(width, height)
# #fig.savefig(path+name+'iapp_bps_kinetics_p1_'+ key+'.png', dpi = 300, bbox_inches='tight')


'''
 Filter the data points with chi^2 less than thr 
'''
  
# thr = 5
# filt = redchi <  thr
# R_filt = R[filt]
# t_filt = t[filt]
# N_filt = N[filt]

# k = 10
# k_filt = 5

# plt.figure()

# fig, ax = plt.subplots()
# fig.subplots_adjust(left=0, bottom=.15, right=.99, top=.97,hspace = 0, wspace = 0)


# y = uniform_filter1d(R, size  = 10)
# y_filt = uniform_filter1d(R_filt, size  = 5)
# plt.plot(t,R,'o', alpha = 0.1, color = 'k', label = r'no $\chi^2$ filter')
# plt.plot(t_filt,R_filt,'o', alpha = 0.2, color = 'b', label = r'$\chi^2 < %.1f$' %(thr))
# plt.plot(t,y, color = 'k', label = 'running average k = %d' %(k))
# plt.plot(t_filt,y_filt, color = 'b', label = r'running average k = %d, $\chi^2 < %.1f$ ' %(k_filt,thr))
# plt.xlabel('time (hrs)')
# plt.ylabel(r'$R_{h}$ (nm)')
# plt.legend()
# print(np.mean(R_filt))

# width = 8
# height = width /1.2

# fig.set_size_inches(width, height)
# fig.savefig('ApoE-Rh-meas-'+ key+'.png', dpi = 300, bbox_inches='tight')



# fig, ax = plt.subplots()
# fig.subplots_adjust(left=0, bottom=.15, right=.99, top=.97,hspace = 0, wspace = 0)

# y = uniform_filter1d(N, size  = 10)
# y_filt = uniform_filter1d(N_filt, size  = 5)
# plt.plot(t,N,'o', alpha = 0.1, color = 'k', label = r'no $\chi^2$ filter')
# plt.plot(t_filt,N_filt,'o', alpha = 0.2, color = 'b', label = r'$\chi^2 < %.1f$' %(thr))
# plt.plot(t,y, color = 'k', label = 'running average k = %d' %(k))
# plt.plot(t_filt,y_filt, color = 'b', label = r'running average k = %d, $\chi^2 < %.1f$ ' %(k_filt,thr))
# plt.xlabel('time (hrs)')
# plt.ylabel(r'$C$ (nM)')
# plt.legend()
# print(np.mean(N_filt))

# width = 8
# height = width /1.2

# fig.set_size_inches(width, height)
# fig.savefig('ApoE-N-meas-'+ key+'.png', dpi = 300, bbox_inches='tight')





# exit()

# plt.figure()
# plt.scatter(R1,redchi1)



# results = []
# ts = np.linspace(2,200,100)
# for t in ts:

#     filt = redchi1 <  t
#     R1_filt = R1[filt]
#     redchi1_filt = redchi1[filt]
#     # plt.scatter(R1_filt,redchi1_filt)

#     # plt.figure()
#     # plt.hist(R1_filt)
#     results.append([np.mean(R1_filt),  np.mean(redchi1_filt), len(R1_filt)])
#     print(np.mean(R1_filt))
    
# results = np.array(results)
    
# plt.figure()
# plt.plot(ts, results[:,0])


