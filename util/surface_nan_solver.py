#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 14:36:46 2019

functions to resolve NaN values near surface

@author: dli84
"""
import numpy as np


def surface_nan_s(data):
    '''
    reslove surface nan for variables that are not u and v
    '''
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0 :
        nan_idx = int(nan_idx[-1] + 1)
        data[:nan_idx] = data[nan_idx+1]
    return(data)

    
def surface_nan_uv(data,z):
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0 :
        nan_idx = int(nan_idx[-1])+1
        para = data[nan_idx]/np.log(z[nan_idx])
        for i in range(0,int(nan_idx)):
            data[i] = para*np.log(z[i])
    return(data)


def surface_nan_w(data):
    '''
    reslove surface nan for vertical veloclity w
    the extra step is taken due to an unkown bug that wrf.interplevel sometimes 
    generates NaN values at certain height
    '''
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0 :
        if nan_idx[-1]> 10 and np.isnan(data[nan_idx[-1]-1])==False and np.isnan(data[nan_idx[-1]+1])==False:
                data[nan_idx[-1]] = data[nan_idx[-1]-1]
        else:
            nan_idx = int(nan_idx[-1] + 1)
            data[:nan_idx] = data[nan_idx+1]
    return(data)