#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# This script contains functions to run create_cfg.py
#
#
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------#


import pandas as pd
from pyproj import Proj
import utm

def domain_location(myProj, centlat, centlon, dx, dy, nx, ny):

    # change lat/lon to UTM to calculate the coordinates of the domain
    utm_centx, utm_centy = myProj(centlon,centlat)
    
    utm_left, utm_right = utm_centx-nx*dx/2, utm_centx+nx*dx/2
    utm_north, utm_south = utm_centy+ny*dy/2, utm_centy-ny*dy/2
    
    # change UTM to lat/lon to locate the domain in the WRF output
    west, north = myProj(utm_left,utm_north,inverse=True)
    east, south = myProj(utm_right,utm_south,inverse=True)
    
    return(west, east, south, north)

def generate_cfg(case_name, dx, dy, dz, nx, ny, nz, west, east, south, north, centlat, centlon):
    cfg = pd.DataFrame({'dx': [dx], 'dy': [dy], 'dz': [dz], 'nx': [nx], 'ny': [ny], 'nz': [nz], 
                        'north': [north], 'south': [south], 'east': [east], 'west': [west], 'centlat':[centlat], 'centlon':[centlon]}) 
    cfg.to_csv('cfg_input/'+ case_name + '.cfg', index=None)
    print('cfg file is ready: '+case_name)
    
def domain_nest(west, south, llx, lly, dx, dy, nx, ny):
    zone = utm.from_latlon(south, west)[2]
    myProj = Proj("+proj=utm +zone=" + str(zone) + ", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    utm_west, utm_south = myProj(west, south)
    utm_cent_ew_d02, utm_cent_sn_d02 = utm_west+llx+(dx*nx/2), utm_south+lly+(dy*ny/2)
    cent_lon_d02, cent_lat_d02 = myProj(utm_cent_ew_d02, utm_cent_sn_d02, inverse=True)
    return(cent_lat_d02, cent_lon_d02)
    
def write_cfg(case_name, dom_dict):
    cfg = pd.DataFrame() 
    for names, values in dom_dict.items():
        cfg[names] = [values]
    cfg.to_csv('cfg_input/'+ case_name + '.cfg', index=None)
    print('cfg file is ready: '+case_name)
    