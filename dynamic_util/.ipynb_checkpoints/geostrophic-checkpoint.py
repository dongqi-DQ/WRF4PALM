#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###------------------------------------------------------------------------###
# WRF4PALM
# Functions to calculate geostrophic wind profiles
# Based on: http://www.meteo.mcgill.ca/~huardda/amelie/geowind.py
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
###------------------------------------------------------------------------###
import numpy as np
def coriolis(lat):  
    """Compute the Coriolis parameter for the given latitude:
    ``f = 2*omega*sin(lat)``, where omega is the angular velocity 
    of the Earth.
    
    Parameters
    ----------
    lat : array
      Latitude [degrees].
    """
    omega   = 7.2921159e-05  # angular velocity of the Earth [rad/s]
    return (2*omega*np.sin(lat/360.0*2*np.pi)) 

def calc_geostrophic_wind(gph, latitude, dy, dx) :
    
    '''
    Use geopotential height to calculte geostrophic wind
    '''
    
    # Set up some constants based on our projection, including the Coriolis parameter and
    # grid spacing, converting lon/lat spacing to Cartesian
    
    f = np.nanmean(coriolis(latitude))
    
    grady = gph[1:, :]- gph[:-1,:]
    gradx = gph[:,:-1] - gph[:, 1:]
    
    
    ug = -np.nanmean(grady/dy * 9.8/f)
    vg = np.nanmean(gradx/dx * 9.8/f)
    

    return(ug, vg)


