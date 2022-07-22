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

def framing_2d_cartesian(lons_wrf,lats_wrf, west,east,south,north, dx_wrf, dy_wrf):
    ## lower left corner
    ll_distance = np.abs(lats_wrf-south) + np.abs(lons_wrf-west)
    nearest_sidx, nearest_widx = np.unravel_index(ll_distance.argmin(),ll_distance.shape)
    if lons_wrf[nearest_sidx, nearest_widx]>west:
        nearest_widx-=1
    if lats_wrf[nearest_sidx, nearest_widx]>south:
        nearest_sidx-=1
        
    ## upper right corner
    ur_distance = np.abs(lats_wrf-north) + np.abs(lons_wrf-east)
    nearest_nidx, nearest_eidx = np.unravel_index(ur_distance.argmin(),ur_distance.shape)
    if lons_wrf[nearest_nidx, nearest_eidx]<east:
        nearest_eidx+=1
    # to make sure the east boundary of PALM domain is within the cropped WRF domain
    if lons_wrf[nearest_nidx, nearest_eidx]-east<dx_wrf:
        nearest_eidx+=1
    if lats_wrf[nearest_nidx, nearest_eidx]<north:
        nearest_nidx+=1
    # to make sure the north boundary of PALM domain is within the cropped WRF domain
    if lats_wrf[nearest_nidx, nearest_eidx]-north<dy_wrf:
        nearest_nidx+=1
    return nearest_widx, nearest_eidx, nearest_sidx, nearest_nidx

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
   
