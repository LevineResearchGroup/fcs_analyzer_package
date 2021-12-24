# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 11:14:36 2021

@author: gwg24
"""

# import numpy as np

# import matplotlib as mpl
# import matplotlib.pyplot as plt
import re
# import lmfit
# from uncertainties import ufloat
# from uncertainties.umath import *

import sys
import platform
sys.path.append("./lib") #point python to location containing the below three modules
# import FCS_fitfunc as ff
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
# import FCS_helpful as fcs

name = 'fluorescin_na_300s_10-28-21'

measurement_group = spt.Read_FCS('/Users/bab226/Documents/yale_research/iapp/fcs/data/BB_Dextran_fcs_control_tests.sptw/' + name)

key = 'DD'

i = 0
for m in measurement_group[:]: #all of them
# for m in measurement_group[:-1]: #exclude the last measurement (e.g., incomplete, or too much evaportation or ...)    
# for m in measurement_group[::5]: #everyother meauserment    

    #SPT64 saves file name in measurement group as e.g., 'ThisName_T1800s_1.ptu'
    #where T1800s means this measurement was taken at 1800 seconds after start.
    #Use regular expression to find pattern Tsomeinteger and then convert it to a number
    measurement_time = float(re.findall(r'T\d+',m['name'])[0][1:]) 
    measurement_time = (measurement_time/60) #now minutes 
    mylabel = 't = %.2f min' %(measurement_time)

    #Use key to decide which correlation curve to look at
    t = m[key]['time']
    G = m[key]['G']
    err = m[key]['err']
    
    t = t/1000
    
    with open('./' + name + '_time_%d.dat' %(measurement_time), "w" ) as f:
        f.write('t,G,err \n')
        for tt,GG,errerr in zip(t, G, err):
            f.write('%e,%e,%e \n' %(tt,GG,errerr))
    
    #Iterator
    i=i+1
    
    # with open('./test%d.dat' %(i), "w" ) as f:
    #     f.write('t \t G \t err \n')
    #     for tt,GG,errerr in zip(t, G, err):
    #         f.write('%e \t %e \t %e \n' %(tt,GG,errerr))
    
    # print(t[0:5])