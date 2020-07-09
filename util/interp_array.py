#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#

# Functions to interpolate array
#
# @initial_author: ricardofaria
# Modified by Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
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
    
    x = np.arange(0, data.shape[0], 1)
    xvals = np.linspace(0, data.shape[0], out_x)
    
    data_res = np.interp(xvals, x, data)

    return (data_res)


