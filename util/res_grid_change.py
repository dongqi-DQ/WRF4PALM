#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:20:12 2017

@author: ricardofaria

change grid resolution
resgrid(data, len(xpoints), len(ypoints))

"""


import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.interpolation import map_coordinates


def grid_points_change(data, out_x, out_y, method) :

    '''

    2d matrix data, x number of points out_x, y number of points out_y, method 'linear' or 'nearest'

    '''
    
#    m = max(data.shape[0], data.shape[1])
#    y = np.linspace(0, 1.0/m, data.shape[0])
#    x = np.linspace(0, 1.0/m, data.shape[1])
#    
#    yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_x), np.linspace(0, 1.0/m, out_y))
    
    y = np.arange(0, data.shape[0], 1)
    x = np.arange(0, data.shape[1], 1)
    interpolating_function = RegularGridInterpolator((y, x), data, method = method)

    yy, xx = np.meshgrid(np.linspace(x[0], x[-1], out_x), np.linspace(y[0], y[-1], out_y))
    
    data_res = interpolating_function((xx, yy))

    return (data_res)


def grid_points_change_1darray(data, out_x) :

    '''

    1d matrix data, x number of points out_x. Output a linear interpolated array

    '''
    
    x = np.arange(0, data.shape[0], 1)
    xvals = np.linspace(0, data.shape[0], out_x)
    
    data_res = np.interp(xvals, x, data)

    return (data_res)


def grid_res_change(data, array_x, array_y, x_res, y_res, method) :

    '''

    2d matrix data, x resolution out_x, y resolution out_y, method 'linear' or 'nearest'

    '''
    
#    array_x = array_x-array_x[0]
#    array_y = array_y-array_y[0]
    
    array_x_res = np.arange(array_x[0], array_x[-1], x_res)
    array_y_res = np.arange(array_y[0], array_y[-1], y_res)
    
#    interpolating_function = RegularGridInterpolator((array_y, array_x), data, method = method)
#    yv, xv = np.meshgrid(array_x_res, array_y_res)
#    data_res = interpolating_function((xv, yv))
    
    X,Y = np.meshgrid(((array_y_res-np.min(array_y))/(np.max(array_y)-np.min(array_y)))*array_y.shape[0],((array_x_res-np.min(array_x))/(np.max(array_x)-np.min(array_x)))*array_x.shape[0])
    #pts = np.asarray(zip(X.ravel(),Y.ravel())).T
    pts = np.array([X.ravel(), Y.ravel()])
    data_res = map_coordinates(data, pts, order = 3, mode = method).reshape(len(array_x_res), len(array_y_res)).T

    return (data_res)
