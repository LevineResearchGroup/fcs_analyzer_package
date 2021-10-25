# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 13:42:18 2021

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

# colors = sns.color_palette("Set2", 8 )
# color1 = colors[0]
# color2 = colors[1]

color2 = 'palevioletred'
color1 = 'lightskyblue'

def label_diff(i,j,text,X,Y):
    x = (X[i]+X[j])/2
    y = 1.05*max(Y[i], Y[j])
    dx = abs(X[i]-X[j])

    props = {'connectionstyle':'bar','arrowstyle':'-',\
                 'shrinkA':20,'shrinkB':20,'linewidth':2}
    ax.annotate(text, xy=(X[i]+dx/2,1.195*y), zorder=10)
    ax.annotate('', xy=(X[i],y), xytext=(X[j],y), arrowprops=props)


def ReadData(filename):
    
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


    return np.array(t), np.array(R), np.array(err)


# def RR_Ratio(R1,E1, R2, E2):
    
#     Z = []
#     EZ = []
#     for r1,e1,r2,e2 in zip(R1,E1,R2,E2):
#         x = ufloat(r1,e1)
#         y = ufloat(r2,e2)
        
#         z = x/y
#         Z.append(z.nominal_value)
#         EZ.append(z.std_dev)

#     return np.array(Z), np.array(EZ)

'''Slow Rh plots'''

Rh_1 = []
Rh_2 = []
width = 3.42
fig = plt.figure(figsize=(width,width/1.62))
ax = fig.add_axes([0, 0, 1, 1])
t, R_DD, err_DD = ReadData('./Results/ADDL_Tube81_Aliquot1_group_Rh2c_DD.dat')
Rh_1.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T81 A1', linestyle = '-', marker = 'o', markersize = 4, color = 'b')

t, R_DD, err_DD = ReadData('./Results/ADDL_Tube81_Aliquot2_Rh2c_DD.dat')
Rh_2.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T81A2', linestyle = '--', marker = 'd', markersize = 4, color = 'b')

t, R_DD, err_DD = ReadData('./Results/ADDL_Tube82_Aliquot1_group_Rh2c_DD.dat')
Rh_1.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T82A1', linestyle = '-', marker = 'o', markersize = 4, color = 'r')
t, R_DD, err_DD = ReadData('./Results/ADDL_Tube82_Aliquot2_Rh2c_DD.dat')
Rh_2.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T82A2', linestyle = '--', marker = 'd', markersize = 4, color = 'r')
t, R_DD, err_DD = ReadData('./Results/ADDL_TubeUnNumb_Aliquot1_group_Rh2c_DD.dat')
Rh_1.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'TxA1', linestyle = '-', marker = 'o', markersize = 4, color = 'k')
t, R_DD, err_DD = ReadData('./Results/ADDL_TubeUnNumb_Aliquot2_Rh2c_DD.dat')
Rh_2.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'TxA2', linestyle = '--', marker = 'd', markersize = 4, color = 'k')

ax.set_ylabel(r'$R_{h}$')
ax.set_xlabel(r't (minutes)')
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.savefig('RhTriplicate', dpi=300, transparent=False, bbox_inches='tight')


Rh_1 = np.array(Rh_1)
Rh_2 = np.array(Rh_2)
labels = ['T81', 'T82', 'TX']
aliquot1_means = Rh_1[:,0]
aliquot1_std = Rh_1[:,1]
aliquot2_means = Rh_2[:,0]
aliquot2_std = Rh_2[:,1]
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

figwidth = 3.42*2
fig = plt.figure(figsize=(figwidth,figwidth/1.62))
# fig, ax = plt.subplots()
ax = fig.add_axes([0, 0, 1, 1])
rects1 = ax.bar(x - width/2, aliquot1_means, yerr = aliquot1_std, width = 0.8*width, label='Aliquot 1', color = color1, alpha = 1, edgecolor = 'k', linewidth = lw)
rects2 = ax.bar(x + width/2, aliquot2_means, yerr = aliquot2_std, width = 0.8*width, label='Aliquot 2', color = color2, alpha = 1, edgecolor = 'k', linewidth = lw)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Rh Slow (nm)')
ax.set_xlabel('Tube')
ax.set_title('Rh slow by Tube and Aliquot')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

ax.bar_label(rects1, padding=3, fmt = '%.2f')
ax.bar_label(rects2, padding=3, fmt = '%.2f')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for tick in ax.get_xticklabels():
    tick.set_rotation(45)

# fig.tight_layout()
plt.savefig('Rhslow_bar.png', dpi=300, transparent=False, bbox_inches='tight')

plt.show()

#is the difference between aliquot 1 and 2 significant?
print('# Rh aliquot differences')
n_a1 = 5
n_a2 = 9
sDelta = np.sqrt( (aliquot1_std**2)/n_a1 + (aliquot2_std**2)/n_a2  ) #Welch's t test, unequal variances, unequal sample sizes
t = (aliquot1_means - aliquot2_means)/sDelta
df = ((sDelta**2)**2) / (  (((aliquot1_std**2)/n_a1)**2)/(n_a1-1) + (((aliquot2_std**2)/n_a2)**2)/(n_a2-1)     )

for tt,ddf in zip(t,df):
    pval = stats.t.sf(np.abs(tt), ddf)*2  # two-sided pvalue = Prob(abs(t)>tt)
    print('t-statistic = %6.3f pvalue = %6.4f' % (tt, pval))


n = 5
print('Aliquot1 Rh')
for i in range(len(aliquot1_means)):
    for j in range(len(aliquot1_means)):
        sDelta = np.sqrt( (aliquot1_std[i]**2)/n + (aliquot1_std[j]**2)/n  ) #Welch's t test, unequal variances, unequal sample sizes
        t = (aliquot1_means[i] - aliquot1_means[j])/sDelta
        df = ((sDelta**2)**2) / (  (((aliquot1_std[i]**2)/n)**2)/(n-1) + (((aliquot1_std[j]**2)/n)**2)/(n-1)     )
        pval = stats.t.sf(np.abs(t), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
        print('%d %d t-statistic = %6.3f pvalue = %6.4f' % (i,j,t, pval))    
        
n = 9
print('Aliquot2 Rh')
for i in range(len(aliquot2_means)):
    for j in range(len(aliquot2_means)):
        sDelta = np.sqrt( (aliquot2_std[i]**2)/n + (aliquot2_std[j]**2)/n  ) #Welch's t test, unequal variances, unequal sample sizes
        t = (aliquot2_means[i] - aliquot2_means[j])/sDelta
        df = ((sDelta**2)**2) / (  (((aliquot2_std[i]**2)/n)**2)/(n-1) + (((aliquot2_std[j]**2)/n)**2)/(n-1)     )
        pval = stats.t.sf(np.abs(t), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
        print('%d %d t-statistic = %6.3f pvalue = %6.4f' % (i,j,t, pval))          







'''
Fraction slow
'''
p1_1 = []
p1_2 = []

width = 3.42
fig = plt.figure(figsize=(width,width/1.62))
ax = fig.add_axes([0, 0, 1, 1])
t, R_DD, err_DD = ReadData('./Results/ADDL_Tube81_Aliquot1_group_p1_2c_DD.dat')
p1_1.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T81A1', linestyle = '-', marker = 'o', markersize = 4, color = 'b')

t, R_DD, err_DD = ReadData('./Results/ADDL_Tube81_Aliquot2_p1_2c_DD.dat')
p1_2.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T81A2', linestyle = '--', marker = 'd', markersize = 4, color = 'b')

t, R_DD, err_DD = ReadData('./Results/ADDL_Tube82_Aliquot1_group_p1_2c_DD.dat')
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T82A1', linestyle = '-', marker = 'o', markersize = 4, color = 'r')
p1_1.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])

t, R_DD, err_DD = ReadData('./Results/ADDL_Tube82_Aliquot2_p1_2c_DD.dat')
p1_2.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T82A2', linestyle = '--', marker = 'd', markersize = 4, color = 'r')

t, R_DD, err_DD = ReadData('./Results/ADDL_TubeUnNumb_Aliquot1_group_p1_2c_DD.dat')
p1_1.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'TxA1', linestyle = '-', marker = 'o', markersize = 4, color = 'k')
t, R_DD, err_DD = ReadData('./Results/ADDL_TubeUnNumb_Aliquot2_p1_2c_DD.dat')
p1_2.append([np.mean(R_DD), np.std(R_DD,ddof = 1)])

ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'TxA2', linestyle = '--', marker = 'd', markersize = 4, color = 'k')
ax.set_ylabel(r'$p_{slow}$')
ax.set_xlabel(r't (minutes)')
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.savefig('FracSlowTriplicate', dpi=300, transparent=False, bbox_inches='tight')

p1_1 = np.array(p1_1)
p1_2 = np.array(p1_2)
labels = ['T81', 'T82', 'TX']
aliquot1_means = p1_1[:,0]
aliquot1_std = p1_1[:,1]
aliquot2_means = p1_2[:,0]
aliquot2_std = p1_2[:,1]
x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

figwidth = 3.42*2
fig = plt.figure(figsize=(figwidth,figwidth/1.62))
# fig, ax = plt.subplots()
ax = fig.add_axes([0, 0, 1, 1])
rects1 = ax.bar(x - width/2, aliquot1_means, yerr = aliquot1_std, width = 0.8*width, label='Aliquot 1', color = color1, alpha = 1, edgecolor = 'k', linewidth = lw)
rects2 = ax.bar(x + width/2, aliquot2_means, yerr = aliquot2_std, width = 0.8*width, label='Aliquot 2', color = color2, alpha = 1, edgecolor = 'k', linewidth = lw)
# label_diff(0,1,'p=0.0370',x-width/2,aliquot1_means)
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Fraction Slow')
ax.set_xlabel('Tube')
ax.set_title('Fraction slow by Tube and Aliquot')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

ax.bar_label(rects1, padding=3, fmt = '%.3f')
ax.bar_label(rects2, padding=3, fmt = '%.3f')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for tick in ax.get_xticklabels():
    tick.set_rotation(45)


# fig.tight_layout()
plt.savefig('fslow_bar.png', dpi=300, transparent=False, bbox_inches='tight')

plt.show()

print('# fslow aliquot differences')
n_a1 = 5
n_a2 = 9
sDelta = np.sqrt( (aliquot1_std**2)/n_a1 + (aliquot2_std**2)/n_a2  ) #Welch's t test, unequal variances, unequal sample sizes
t = (aliquot1_means - aliquot2_means)/sDelta
df = ((sDelta**2)**2) / (  (((aliquot1_std**2)/n_a1)**2)/(n_a1-1) + (((aliquot2_std**2)/n_a2)**2)/(n_a2-1)     )

for tt,ddf in zip(t,df):
    pval = stats.t.sf(np.abs(tt), ddf)*2  # two-sided pvalue = Prob(abs(t)>tt)
    print('t-statistic = %6.3f pvalue = %6.4f' % (tt, pval))

#within aliquot difference b/w tubes
#aliquot 1

aliquot1_means_joined = np.hstack((aliquot1_means, aliquot1_means[0]))
print(aliquot1_means_joined)

n = 5
print('Aliquot1 fslow')
for i in range(len(aliquot1_means)):
    for j in range(len(aliquot1_means)):
        sDelta = np.sqrt( (aliquot1_std[i]**2)/n + (aliquot1_std[j]**2)/n  ) #Welch's t test, unequal variances, unequal sample sizes
        t = (aliquot1_means[i] - aliquot1_means[j])/sDelta
        df = ((sDelta**2)**2) / (  (((aliquot1_std[i]**2)/n)**2)/(n-1) + (((aliquot1_std[j]**2)/n)**2)/(n-1)     )
        pval = stats.t.sf(np.abs(t), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
        print('%d %d t-statistic = %6.3f pvalue = %6.4f' % (i,j,t, pval))    
        
n = 9
print('Aliquot2 fslow')
for i in range(len(aliquot2_means)):
    for j in range(len(aliquot2_means)):
        sDelta = np.sqrt( (aliquot2_std[i]**2)/n + (aliquot2_std[j]**2)/n  ) #Welch's t test, unequal variances, unequal sample sizes
        t = (aliquot2_means[i] - aliquot2_means[j])/sDelta
        df = ((sDelta**2)**2) / (  (((aliquot2_std[i]**2)/n)**2)/(n-1) + (((aliquot2_std[j]**2)/n)**2)/(n-1)     )
        pval = stats.t.sf(np.abs(t), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
        print('%d %d t-statistic = %6.3f pvalue = %6.4f' % (i,j,t, pval))          

'''
N
'''


width = 3.42
fig = plt.figure(figsize=(width,width/1.62))
ax = fig.add_axes([0, 0, 1, 1])
t, R_DD, err_DD = ReadData('./Results/ADDL_Tube81_Aliquot1_group_N_2c_DD.dat')
print('Mean N = %.4f +/- %.5f' %(np.mean(R_DD), np.std(R_DD, ddof =1)))
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T81A1', linestyle = '-', marker = 'o', markersize = 4, color = 'b')
t, R_DD, err_DD = ReadData('./Results/ADDL_Tube81_Aliquot2_N_2c_DD.dat')
print('Mean N = %.4f +/- %.5f' %(np.mean(R_DD), np.std(R_DD, ddof =1)))
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T81A2', linestyle = '--', marker = 'd', markersize = 4, color = 'b')

t, R_DD, err_DD = ReadData('./Results/ADDL_Tube82_Aliquot1_group_N_2c_DD.dat')
print('Mean N = %.4f +/- %.5f' %(np.mean(R_DD), np.std(R_DD, ddof =1)))
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T82A1', linestyle = '-', marker = 'o', markersize = 4, color = 'r')
t, R_DD, err_DD = ReadData('./Results/ADDL_Tube82_Aliquot2_N_2c_DD.dat')
print('Mean N = %.4f +/- %.5f' %(np.mean(R_DD), np.std(R_DD, ddof =1)))
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'T82A2', linestyle = '--', marker = 'd', markersize = 4, color = 'r')


t, R_DD, err_DD = ReadData('./Results/ADDL_TubeUnNumb_Aliquot1_group_N_2c_DD.dat')
print('Mean N = %.4f +/- %.5f' %(np.mean(R_DD), np.std(R_DD, ddof =1)))
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'TxA1', linestyle = '-', marker = 'o', markersize = 4, color = 'k')
t, R_DD, err_DD = ReadData('./Results/ADDL_TubeUnNumb_Aliquot2_N_2c_DD.dat')
print('Mean N = %.4f +/- %.5f' %(np.mean(R_DD), np.std(R_DD, ddof =1)))
ax.errorbar(t, R_DD, yerr = err_DD, linewidth =1, label = 'TxA2', linestyle = '--', marker = 'd', markersize = 4, color = 'k')

ax.set_ylabel(r'$N$')
ax.set_xlabel(r't (minutes)')
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.savefig('NTriplicate', dpi=300, transparent=False, bbox_inches='tight')