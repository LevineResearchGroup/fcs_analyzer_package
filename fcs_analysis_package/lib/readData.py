#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 18:40:35 2021

Written by: gwg24

Purpose: Read functions for secondary analysis of Batch_Fitting results.

"""
#Packages to import
import csv
import numpy as np

#Read Data from Batch_Fitting_FCS python script.
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

