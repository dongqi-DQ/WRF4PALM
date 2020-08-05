#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# WRF4PALM
# This script generates a cfg input file for create_dynamic.py
# This file requires users provide:
#     latitude and longitude of domain centre
#     domain size and resolution (dx, dy, dz, nx, ny, nz)
#
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------#

import numpy as np
import pandas as pd
from pyproj import Proj
import utm

from util.loc_dom import *

###############################################################################
##----------------------------- Domain Config -------------------------------##
###############################################################################

# domain 01
case_name_d01 = 'chch_NW_10m'
centlat_d01 =-43.487
centlon_d01 = 172.537
dx_d01 = 10
dy_d01 = 10
dz_d01 = 10
nx_d01 = 360
ny_d01 = 360
nz_d01 = 200

# calculate zone of projection for lat/lon calculation
zone = utm.from_latlon(centlat_d01,centlon_d01)[2]

myProj = Proj("+proj=utm +zone=" + str(zone) + ", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

west_d01, east_d01, south_d01, north_d01 = domain_location(myProj, centlat_d01, centlon_d01, dx_d01, dy_d01, nx_d01, ny_d01)
# write a cfg file
generate_cfg(case_name_d01, dx_d01, dy_d01, dz_d01, nx_d01, ny_d01, nz_d01, west_d01, east_d01, south_d01, north_d01, centlat_d01, centlon_d01)

###############################################################################
##----------------------------- Nested Domains ------------------------------##
###############################################################################
## domain 02 (This is only a reference for nested domains)
case_name_d02 = 'chch_NW_5m'

# lower left x and y in meters
ll_x_d02 = 1500
ll_y_d02 = 1500

dx_d02 = 5
dy_d02 = 5
dz_d02 = 5
nx_d02 = 120
ny_d02 = 120
nz_d02 = 120

centlat_d02, centlon_d02 = domain_nest(west_d01, south_d01, ll_x_d02, ll_y_d02, dx_d02, dy_d02, nx_d02, ny_d02)

west_d02, east_d02, south_d02, north_d02 = domain_location(myProj, centlat_d02, centlon_d02, dx_d02, dy_d02, nx_d02, ny_d02)

generate_cfg(case_name_d02, dx_d02, dy_d02, dz_d02, nx_d02, ny_d02, nz_d02, west_d02, east_d02, south_d02, north_d02, centlat_d02, centlon_d02)


