#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 Created on Tue Nov 19 11:17:36 2019
 Python script to generate a cfg file for WRF-PALM coupler to locate the interpolation
 domain. 
 This script requires users provide:
    latitude and longitude of domain centre
    domain size and resolution (dx, dy, dz, nx, ny, nz)

 
@author: Dongqi Lin
 To do list:
     add vertical streching feature
'''


import pandas as pd
from pyproj import Proj
import utm

def domain_locating(midlat, midlon, dx, dy, nx, ny):
    
    # change lat/lon to UTM to calculate the coordinates of the domain
    utm_midx, utm_midy = myProj(midlon,midlat)
    
    utm_left, utm_right = utm_midx-nx*dx/2, utm_midx+nx*dx/2
    utm_north, utm_south = utm_midy+ny*dy/2, utm_midy-ny*dy/2
    
    # change UTM to lat/lon to locate the domain in the WRF output
    west, north = myProj(utm_left,utm_north,inverse=True)
    east, south = myProj(utm_right,utm_south,inverse=True)
    
    return(west, east, south, north)

def generate_cfg(case_name, dx, dy, dz, nx, ny, nz, west, east, south, north, midlat, midlon):
    cfg = pd.DataFrame({'dx': [dx], 'dy': [dy], 'dz': [dz], 'nx': [nx], 'ny': [ny], 'nz': [nz], 
                        'north': [north], 'south': [south], 'east': [east], 'west': [west], 'midlat':[midlat], 'midlon':[midlon]}) 
    cfg.to_csv('cfg_input/'+ case_name + '.cfg', index=None)
    print('cfg file is ready: '+case_name)
    
def dom_lat_lon(west, south, llx, lly, dx, dy, nx, ny):
    zone = utm.from_latlon(south, west)[2]
    myProj = Proj("+proj=utm +zone=" + str(zone) + ", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    utm_west, utm_south = myProj(west, south)
    utm_mid_ew_d02, utm_mid_sn_d02 = utm_west+llx+(dx*nx/2), utm_south+lly+(dy*ny/2)
    mid_lon_d02, mid_lat_d02 = myProj(utm_mid_ew_d02, utm_mid_sn_d02, inverse=True)
    return(mid_lat_d02, mid_lon_d02)

###############################################################################
##----------------------------- Define domains ------------------------------##
###############################################################################

# domain 01
case_name_d01 = 'chch_50m_example'
midlat_d01 =-43.508760 
midlon_d01 = 172.664099
dx_d01 = 50
dy_d01 = 50
dz_d01 = 50
nx_d01 = 96
ny_d01 = 96
nz_d01 = 96


zone = utm.from_latlon(midlat_d01,midlon_d01)[2]#'32'  # this is the UTM zone for projection 
#32N for Netherland; 33U for Berlin; 59S for CHCH

myProj = Proj("+proj=utm +zone=" + str(zone) + ", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

west_d01, east_d01, south_d01, north_d01 = domain_locating(midlat_d01, midlon_d01, dx_d01, dy_d01, nx_d01, ny_d01)

generate_cfg(case_name_d01, dx_d01, dy_d01, dz_d01, nx_d01, ny_d01, nz_d01, west_d01, east_d01, south_d01, north_d01, midlat_d01, midlon_d01)


'''
The code below allows users to enable PALM self nesting and calculate the coordinates information
for the nested domain. The offiline nesting is not needed for the nested domains. I keep this feature just 
to keep a reference of the nested domains.
'''
### domain 02
#case_name_d02 = 'chch_50m_example_d02'
#
## lower left x and y in meters
#ll_x_d02 = 1000
#ll_y_d02 = 1000
#
#dx_d02 = 10
#dy_d02 = 10
#dz_d02 = 10
#nx_d02 = 96
#ny_d02 = 96
#nz_d02 = 96
#
#midlat_d02, midlon_d02 = dom_lat_lon(west_d01, south_d01, ll_x_d02, ll_y_d02, dx_d02, dy_d02, nx_d02, ny_d02)
#
#west_d02, east_d02, south_d02, north_d02 = domain_locating(midlat_d02, midlon_d02, dx_d02, dy_d02, nx_d02, ny_d02)
#
#generate_cfg(case_name_d02, dx_d02, dy_d02, dz_d02, nx_d02, ny_d02, nz_d02, west_d02, east_d02, south_d02, north_d02, midlat_d02, midlon_d02)


