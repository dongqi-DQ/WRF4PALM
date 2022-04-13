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
import numpy as np

def framing_2d_cartesian(lons_wrf,lats_wrf, west,east,south,north):
    nearest_wix = np.argmin(np.abs(lons_wrf[0]-west))
    if lons_wrf[0][nearest_wix]>west:
        nearest_wix-=1
    nearest_eix = np.argmin(np.abs(lons_wrf[0]-east))
    if lons_wrf[0][nearest_eix]<east:
        nearest_eix+=1

    lats_wrf_y = np.array([x[0] for x in lats_wrf])

    nearest_six = np.argmin(np.abs(lats_wrf_y-south))
    if lats_wrf_y[nearest_six]>south:
        nearest_six-=1
    nearest_nix = np.argmin(np.abs(lats_wrf_y-north))
    if lats_wrf_y[nearest_nix]<north:
        nearest_nix+=1

    return nearest_wix, nearest_eix, nearest_six, nearest_nix


def nearest_2d(array, number):
    
    '''
    
    find nearest index value and index in a 2D array.
    
    nearest(array, number)
    
    return(nearest_number, nearest_index)
    
    '''
      
    index = np.where(np.abs(array-number) == np.nanmin(np.abs(array-number)))#index is single pair of naturals, if only array isn't anywhere cartesian.
    idx, idy = index[0], index[1]
    nearest_number = array[idx, idy]#1D array of length 1 or Nx of Ny.
    nearest_index = [idx[0], idy[0]]
    
    return(nearest_number[0], nearest_index)

def nearest_1d(array, number):
    
    '''
    
    find nearest index value and index in array.
    
    nearest(array, number)
    
    return(nearest_number, nearest_index)
    
    '''
    
    nearest_index = np.where(np.abs(array-number) == np.nanmin(np.abs(array-number)))
    nearest_index = int(nearest_index[0][0])
    nearest_number = array[nearest_index]
    
    return(nearest_number, nearest_index)
   
