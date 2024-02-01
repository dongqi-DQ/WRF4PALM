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
import numpy as np
from pyproj import Proj,  Transformer, CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info
from math import ceil, floor
import numpy as np

def domain_location(palm_proj, wgs_proj, centlat, centlon, dx, dy, nx, ny):
    '''
    Identify west, east, south, and north boundaries of doamin based on lat/lon at
    domain centre
    '''
    if palm_proj == wgs_proj:
        # convert to UTM to calculate the lat/lon
        utm_crs_list = query_utm_crs_info(
            datum_name="WGS 84",
            area_of_interest=AreaOfInterest(
                west_lon_degree=centlon,
                south_lat_degree=centlat,
                east_lon_degree=centlon,
                north_lat_degree=centlat,
            ),
        )
        utm_crs = CRS.from_epsg(utm_crs_list[0].code)
        trans_utm2wgs = Transformer.from_proj( utm_crs, wgs_proj, always_xy=True)
        trans_wgs2utm = Transformer.from_proj( wgs_proj, utm_crs, always_xy=True)

        centx, centy =  trans_wgs2utm.transform(centlon, centlat)
        west_utm, east_utm   = centx-nx*dx/2, centx+nx*dx/2
        north_utm, south_utm = centy+ny*dy/2, centy-ny*dy/2
        west, north = trans_utm2wgs.transform( west_utm,north_utm)
        east, south = trans_utm2wgs.transform( east_utm,south_utm)

    else:
    # change lat/lon to UTM to calculate the coordinates of the domain

        inProj =  wgs_proj
        outProj = palm_proj
        transformer = Transformer.from_proj(inProj, outProj, always_xy=True)
        centx, centy = transformer.transform(centlon,centlat)

        west, east   = centx-nx*dx/2, centx+nx*dx/2
        north, south = centy+ny*dy/2, centy-ny*dy/2

    return west, east, south, north, centx, centy


def generate_cfg(case_name, dx, dy, dz, nx, ny, nz, west, east, south, north, centlat, centlon, z_origin):
    cfg = pd.DataFrame({'dx': [dx], 'dy': [dy], 'dz': [dz], 'nx': [nx], 'ny': [ny], 'nz': [nz], 'z_origin': [z_origin],
                        'north': [north], 'south': [south], 'east': [east], 'west': [west], 'centlat':[centlat], 'centlon':[centlon]})
    cfg.to_csv('cfg_files/'+ case_name + '.cfg', index=None)
    print('cfg file is saved: '+case_name)



def calc_stretch(z, dz, zw, dz_stretch_factor, dz_stretch_level, dz_max):
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
