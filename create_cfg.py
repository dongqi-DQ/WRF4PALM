#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# This script is the cfg file generator for WRF-PALM coupler
# This file requires users provide:
#     latitude and longitude of domain centre
#     domain size and resolution (dx, dy, dz, nx, ny, nz)
#
# To do list:
#     add vertical streching feature
#
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------#

import numpy as np
import pandas as pd
from pyproj import Proj
import utm

from dynamic_util.loc_dom import *

###############################################################################
##----------------------------- Define domains ------------------------------##
###############################################################################

# domain 01
case_name_d01 = 'chch_airport_256m'
centlat_d01 =-43.493 
centlon_d01 = 172.537
dx_d01 = 256
dy_d01 = 256
dz_d01 = 256
nx_d01 = 80
ny_d01 = 64
nz_d01 = 64


zone = utm.from_latlon(centlat_d01,centlon_d01)[2]

myProj = Proj("+proj=utm +zone=" + str(zone) + ", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

west_d01, east_d01, south_d01, north_d01 = domain_location(myProj, centlat_d01, centlon_d01, dx_d01, dy_d01, nx_d01, ny_d01)

generate_cfg(case_name_d01, dx_d01, dy_d01, dz_d01, nx_d01, ny_d01, nz_d01, west_d01, east_d01, south_d01, north_d01, centlat_d01, centlon_d01)


## domain 02 (a reference for nested domains)
case_name_d02 = 'chch_airport_128m'

# lower left x and y in meters
ll_x_d02 = 2560
ll_y_d02 = 4096

dx_d02 = 128
dy_d02 = 128
dz_d02 = 128
nx_d02 = 120
ny_d02 = 64
nz_d02 = 96

centlat_d02, centlon_d02 = domain_nest(west_d01, south_d01, ll_x_d02, ll_y_d02, dx_d02, dy_d02, nx_d02, ny_d02)

west_d02, east_d02, south_d02, north_d02 = domain_location(myProj, centlat_d02, centlon_d02, dx_d02, dy_d02, nx_d02, ny_d02)

generate_cfg(case_name_d02, dx_d02, dy_d02, dz_d02, nx_d02, ny_d02, nz_d02, west_d02, east_d02, south_d02, north_d02, centlat_d02, centlon_d02)


## domain 03 (a reference for nested domains)
case_name_d03 = 'chch_airport_64m'
#
## lower left x and y in meters
#ll_x_d03 = 5120
#ll_y_d03 = 5120
#
#dx_d03 = 64
#dy_d03 = 64
#dz_d03 = 64
#nx_d03 = 160
#ny_d03 = 96
#nz_d03 = 128
#
#centlat_d03, centlon_d03 = domain_nest(west_d01, south_d01, ll_x_d03, ll_y_d03, dx_d03, dy_d03, nx_d03, ny_d03)
#
#west_d03, east_d03, south_d03, north_d03 = domain_location(myProj, centlat_d03, centlon_d03, dx_d03, dy_d03, nx_d03, ny_d03)
#
#generate_cfg(case_name_d03, dx_d03, dy_d03, dz_d03, nx_d03, ny_d03, nz_d03, west_d03, east_d03, south_d03, north_d03, centlat_d03, centlon_d03)

