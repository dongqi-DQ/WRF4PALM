#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# WRF4PALM
# Functions to interpolate 2d and 1d array
# (adopted from WRF2PALM)
# @ WRF2PALM author: Ricardo Faria (https://github.com/ricardo88faria/WRF2PALM)
# @ WRF4PALM contact: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------#



import numpy as np
from scipy.interpolate import RegularGridInterpolator


def interp_array_2d(data, out_x, out_y, method) :

    '''

    2d matrix data, x number of points out_x, y number of points out_y, method 'linear' or 'nearest'

    '''
    
    y = np.arange(0, data.shape[0], 1)
    x = np.arange(0, data.shape[1], 1)
    interpolating_function = RegularGridInterpolator((y, x), data, method = method)

    yy, xx = np.meshgrid(np.linspace(x[0], x[-1], out_x), np.linspace(y[0], y[-1], out_y))
    
    data_res = interpolating_function((xx, yy))

    return (data_res)


def interp_array_1d(data, out_x) :

    '''

    1d matrix data, x number of points out_x. Output a linear interpolated array

    '''
    
def interp_array_1d(data, out_x) :

    '''

    1d matrix data, x number of points out_x. Output a linear interpolated array

    '''
    
    x = np.linspace(0, 1, data.shape[0])
    xvals = np.linspace(0, 1, out_x)
    
    data_res = np.interp(xvals, x, data)

    return (data_res)


