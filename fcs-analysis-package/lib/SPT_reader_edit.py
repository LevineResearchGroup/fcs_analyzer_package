# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 15:01:27 2021

@author: gwg24

w/ catch for bug in SPT output of xcor
"""
import csv
import numpy as np

# import matplotlib as mpl
# import matplotlib.pyplot as plt
import re
# import lmfit


def Read_FCS(filename):
    
    #Convert encoding from PC for linux/Mac.
    coding1 = "iso-8859-1"
    coding2 = "utf-8"

    f= open(filename + ".dat", 'r', encoding=coding1)
    content= f.read()
    f.close()
    
    filename_convert = filename[:-4]+'.tmp'
    f= open(filename_convert, 'w', encoding=coding2)
    f.write(content)
    f.close()

    print("done")

    def multipliercheck(col):
        #Annoying feature of SPT export is that sometimes G or err will be exported as 10^-3 and sometimes not.
        #this will check header for pattern "minus sign followed by integer")
        test_str = re.findall(r'-\d+',header[1][col])
        if len(test_str) > 0:
            M = 10**(int(test_str[0]))
        else:
            M = 1 
        # print(M) 
        return M
    
    hl = 2 #number of header line
    header = [] #an empty list to which we will append header lines

    #First read file to get header information, and number of measurements
    with open(filename_convert, newline = '') as data_file:                                                                                          
        reader = csv.reader(data_file, delimiter='\t')
        
        for i in range(hl):
            #Read and store headerlines (non-data)
            header.append(next(reader))


    if header[1][2][:4] == '±Err':    
        #number of measurements in group is number of rows/9 (two autocorrelations, one cross correlation and each correlation has  time, G, err, i.e., 3*3 = 9)
        ncurves = len(header[0]) 
        n_measurement= ncurves/9
        # assert n_measurement.is_integer()
        print(ncurves)
        print('Xcor format identified as t,G, err')    
        print('Reading .dat files comprised of %d measurements' %(n_measurement))
        xcor_flag = False
    else:
        #number of measurements in group is number of rows/9 (two autocorrelations, one cross correlation and each correlation has  time, G, err, i.e., 3*3 = 9)
        ncurves = len(header[0]) 
        print(ncurves)
        n_measurement= (ncurves+1)/9
        # assert n_measurement.is_integer()
        print('Xcor format identified as t,G, no error')        
        print('Reading .dat files comprised of %d measurements' %(n_measurement))
        xcor_flag = True

    #Populate with data by iteratively reading file
    measurement_group = []
    print('Measurement names are: ....')
    for i in range(int(n_measurement)):        
        measurement = {} #initialize a dictionary to store measurement


        #0, 9, 18
        if xcor_flag:
            ref_col = 8*i
        else:
            ref_col = 9*i #time column will be a reference colomn, from whihc other column positions are deduced.

        measurement['name'] = header[0][0 + ref_col]
        print(measurement['name'])
        
        tx = [] 
        Gx = [] 
        errx = []
        tdd = [] 
        Gdd = [] 
        errdd = []
        taa = [] 
        Gaa = [] 
        erraa = []

        with open(filename_convert, newline = '') as data_file:                                                                                          
            reader = csv.reader(data_file, delimiter='\t')
            for j in range(hl):
                #just skip headerlines, already stored
                next(reader)
            for row in reader:
                #Read remaining lines (row at a time)
                
                tx.append(float(row[0 + ref_col]))
                Gx.append(float(row[1 + ref_col]))
                
                if xcor_flag:
                    errx.append(float(1)) #1 so weights = 1/err = 1
                    tdd.append(float(row[2 + ref_col]))
                    Gdd.append(float(row[3 + ref_col]))
                    errdd.append(float(row[4 + ref_col]))
                    
                    taa.append(float(row[5 + ref_col]))
                    Gaa.append(float(row[6 + ref_col]))
                    erraa.append(float(row[7 + ref_col]))
                else:
                    errx.append(float(row[2 + ref_col]))
                
                    tdd.append(float(row[3 + ref_col]))
                    Gdd.append(float(row[4 + ref_col]))
                    errdd.append(float(row[5 + ref_col]))
                    
                    taa.append(float(row[6 + ref_col]))
                    Gaa.append(float(row[7 + ref_col]))
                    erraa.append(float(row[8 + ref_col]))
            
        
        curve_x = {}
        curve_x['time'] = np.asarray(tx)
        curve_x['G'] = np.asarray(Gx)*multipliercheck(1 + ref_col)
        curve_x['err'] = np.asarray(errx)*multipliercheck(2 + ref_col)
        measurement['DxA'] = curve_x
        
        curve_dd = {}
        curve_dd['time'] = np.asarray(tdd)
        curve_dd['G'] = np.asarray(Gdd)*multipliercheck(4 +ref_col)
        curve_dd['err'] = np.asarray(errdd)*multipliercheck(5 + ref_col)
        measurement['DD'] = curve_dd
        
        curve_aa = {}
        curve_aa['time'] = np.asarray(taa)
        curve_aa['G'] = np.asarray(Gaa)*multipliercheck(7 + ref_col)
        curve_aa['err'] = np.asarray(erraa)*multipliercheck(8 + ref_col)
        measurement['AA'] = curve_aa
        
        
        measurement_group.append(measurement)
    
    return measurement_group