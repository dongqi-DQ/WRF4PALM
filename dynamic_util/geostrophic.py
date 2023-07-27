#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###------------------------------------------------------------------------###
# WRF4PALM
# Functions to calculate geostrophic wind profiles
# Based on: http://www.meteo.mcgill.ca/~huardda/amelie/geowind.py
# @author: Dongqi Lin (dongqi.lin@canterbury.ac.nz)
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
    import numpy as np
    omega   = 7.2921159e-05  # angular velocity of the Earth [rad/s]
    return (2*omega*np.sin(lat/360.0*2*np.pi)) 

def rho(T, p):
    
    """
    
    Calculates air density (rho)
    
    """
    
    
    Rd = 287.0

#    Tv   = T * (1+0.61*qv) # Virtual temperature

    rho = p / (Rd * T) # Air density [kg m^-3]

    return(rho)
    
def __midpoints_1d(a):
    """Return `a` linearly interpolated at the mid-points."""
    return((a[:-1] + a[1:])/2.0)
    
def midpoints(a,  axis=None):
    """Return `a` linearly interpolated at the mid-points.
    
    Parameters
    ----------
    a : array-like 
      Input array.
    axis : int or None
      Axis along which the interpolation takes place. None stands for all axes. 
    
    Returns
    -------
    out : ndarray 
      Input array interpolated at the midpoints along the given axis. 
      
    Examples
    --------
    >>> a = [1,2,3,4]
    >>> midpoints(a)
    array([1.5, 2.5, 3.5])
    """
    import numpy as np
    x = np.asarray(a)
    if axis is not None:
        return(np.apply_along_axis(__midpoints_1d,  axis, x))
    else:
        for i in range(x.ndim):
            x = midpoints(x,  i)
        return(x)
    
def calc_geostrophic_wind_plevels(array_2d_press, array_2d_temp, array_1d_lat, array_1d_lon,dy, dx) :
    
    '''

    Calculate Geostrophic wind profile (1 point value representing input 2d array area).
    Based on Practical_Meteorology-v1.02b-WholeBookColor pag.302
    
    Parameters
    ----------
    array_2d_press : read numpy array [Pa]
    array_2d_temp : read numpy array [k]
    array_1d_lat : read numpy array [deg]
    array_1d_lon : read numpy array [deg]

    Returns
    -------
    array : return interplated and extrapolated value

    '''
  
    
    # Set up some constants based on our projection, including the Coriolis parameter and
    # grid spacing, converting lon/lat spacing to Cartesian
    
    fx = np.nanmean(coriolis(array_1d_lat))*np.mean(dx)
    fy = np.nanmean(coriolis(array_1d_lat))*np.mean(dy)
    
    
    rho_tmp = np.nanmean(rho(array_2d_temp, array_2d_press))
       
    
    gradx = np.zeros((len(array_1d_lat),len(array_1d_lon)))
    grady = np.zeros((len(array_1d_lat),len(array_1d_lon)))
    
    for i in range(0,len(array_1d_lon)-1):
        gradx[:,i] = array_2d_press[:,i+1]-array_2d_press[:,i] 
    
    for j in range(0,len(array_1d_lat)-1):
        grady[j,:] = array_2d_press[j+1,:]-array_2d_press[j,:] 
 
    gradx = midpoints(gradx,  axis=0)
    grady = midpoints(grady,  axis=1)
    
    ug_tmp = np.nanmean((-1/(rho_tmp*fy))*grady)
    vg_tmp = np.nanmean((1/(rho_tmp*fx))*gradx)
    
    
    geo_wind = np.array([ug_tmp, vg_tmp])
    

    return(geo_wind)



def calc_geostrophic_wind_zlevels(gph, latitude, dy, dx) :
    
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