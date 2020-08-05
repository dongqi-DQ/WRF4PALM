#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# WRF4PALM
#
# function to find nearest index value and index in array
# (adopted from WRF2PALM)
# @ WRF2PALM author: Ricardo Faria (https://github.com/ricardo88faria/WRF2PALM)
# @ WRF4PALM contact: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------#


def nearest(array, number):
    
    '''
    
    find nearest index value and index in array.
    
    nearest(array, number) 
    
    return(nearest_number, nearest_index)
    
    '''
    
    import numpy as np
    
    nearest_index = np.where(np.abs(array-number) == np.nanmin(np.abs(array-number)))
    nearest_index = int(nearest_index[0])
    nearest_number = array[nearest_index]
    
    return(nearest_number, nearest_index)
    