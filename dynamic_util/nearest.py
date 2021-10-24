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
def nearest_2d(array, number):
    
    '''
    
    find nearest index value and index in a 2D array.
    
    nearest(array, number) 
    
    return(nearest_number, nearest_index)
    
    '''
    
    import numpy as np
    
    index = np.where(np.abs(array-number) == np.nanmin(np.abs(array-number)))
    idx, idy = index[0], index[1]
    nearest_number = array[idx, idy]
    nearest_index = [idx[0], idy[0]]
    
    return(nearest_number[0], nearest_index)

def nearest(array, number):
    
    '''
    
    find nearest index value and index in array.
    
    nearest(array, number) 
    
    return(nearest_number, nearest_index)
    
    '''
    
    import numpy as np
    
    nearest_index = np.where(np.abs(array-number) == np.nanmin(np.abs(array-number)))
    nearest_index = int(nearest_index[0][0])
    nearest_number = array[nearest_index]
    
    return(nearest_number, nearest_index)
   
