#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# WRF4PALM
# This script contains functions for create_cfg.py
#   - function to calculate lat/lon
#   - function to write cfg file
#   - function to calculate centre lat/lon for nested domains
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------#

import pandas as pd
from pyproj import Proj,  Transformer

def domain_location(palm_proj, wgs_proj, centlat, centlon, dx, dy, nx, ny):
    '''
    Identify west, east, south, and north boundaries of doamin based on lat/lon at
    domain centre 
    '''
    if palm_proj == wgs_proj:
        # no conversion needed when latlong projection is used
        centx =  centlon
        centy =  centlat
    else:
    # change lat/lon to UTM to calculate the coordinates of the domain

        inProj =  wgs_proj
        outProj = palm_proj
        transformer = Transformer.from_proj(inProj, outProj)
        centx, centy = transformer.transform(centlon,centlat)
    
    west, east   = centx-nx*dx/2, centx+nx*dx/2
    north, south = centy+ny*dy/2, centy-ny*dy/2
    
    return west, east, south, north


def generate_cfg(case_name, dx, dy, dz, nx, ny, nz, west, east, south, north, centlat, centlon, z_origin):
    cfg = pd.DataFrame({'dx': [dx], 'dy': [dy], 'dz': [dz], 'nx': [nx], 'ny': [ny], 'nz': [nz], 'z_origin': [z_origin],
                        'north': [north], 'south': [south], 'east': [east], 'west': [west], 'centlat':[centlat], 'centlon':[centlon]}) 
    cfg.to_csv('cfg_files/'+ case_name + '.cfg', index=None)
    print('cfg file is saved: '+case_name)
    
    

def calc_stretch(z, dz):
    '''
    vertically streched levels
    '''
    dz_lvl = np.zeros_like(z)
    dz_lvl[:] = dz
    z_stretch = np.copy(z)
    zw_stretch = np.copy(zw)
    for idz, height in enumerate(z):
        if height>dz_stretch_level:
            dz_lvl[idz] = dz_lvl[idz-1]*dz_stretch_factor
            if dz_lvl[idz]<=dz_max:
                z_stretch[idz] = z_stretch[idz-1]+dz_lvl[idz]
            else:
                z_stretch[idz] = z_stretch[idz-1]+dz_max
    for i in range(0, zw.shape[0]):
        zw_stretch[i] = (z_stretch[i]+z_stretch[i+1])*0.5
    return(z_stretch, zw_stretch)