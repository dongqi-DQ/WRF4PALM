#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:29:13 2018

@author: ricardofaria

find nearest index value and index in array

"""

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
    