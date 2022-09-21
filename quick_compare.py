#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------
# WRF4PALM quick compare
#--------------------------------------------------------------------------------
# quickly compare WRF output with the dynamic driver
# this script offers 3 types of comparison
# - zcross: vertical cross sections at west/east/south/north boundaries of PALM domain
#           for the given variable and timestamp
# - pr: vertical profiles (horizontal mean is taken) at west/east/south/north
#       boundaries of PALM domain for the given variable and timestamp
# - ts: time series horizontal mean is taken) at west/east/south/north
#       boundaries of PALM domain for the given variable and altitude
# How to use:
# 1. python quick_compare.py namelist.wrf4palm [plot type] variable
# 2. the script will ask for the required timestamp or altitude depending on the plot type
#
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------

import sys
import salem
import xarray as xr
from pyproj import Proj, Transformer
import configparser
import ast
from glob import glob
import numpy as np
from math import ceil, floor
from datetime import datetime, timedelta
from dynamic_util.nearest import framing_2d_cartesian, nearest_1d
from dynamic_util.loc_dom import calc_stretch, domain_location, generate_cfg
from dynamic_util.surface_nan_solver import *
import warnings
## supress warnings
## switch to other actions if needed
warnings.filterwarnings("ignore", '.*pyproj.*')

import matplotlib.pyplot as plt
import pandas as pd

params={
    'axes.labelsize'  : '16',
    'axes.titlesize'  : '16',
    'xtick.labelsize' :'14',
    'ytick.labelsize' :'16',
    'lines.linewidth' : '2' ,
    'legend.fontsize' : '16',
    'figure.figsize'   : '16, 9'
}
plt.rcParams.update(params)


settings_cfg = configparser.ConfigParser(inline_comment_prefixes='#')
config = configparser.RawConfigParser()
config.read(sys.argv[1])# read namelist
case_name =  ast.literal_eval(config.get("case", "case_name"))[0]
max_pool =  ast.literal_eval(config.get("case", "max_pool"))[0]

palm_proj_code = ast.literal_eval(config.get("domain", "palm_proj"))[0]
centlat = ast.literal_eval(config.get("domain", "centlat"))[0]
centlon = ast.literal_eval(config.get("domain", "centlon"))[0]
dx = ast.literal_eval(config.get("domain", "dx"))[0]
dy = ast.literal_eval(config.get("domain", "dy"))[0]
dz = ast.literal_eval(config.get("domain", "dz"))[0]
nx = ast.literal_eval(config.get("domain", "nx"))[0]
ny = ast.literal_eval(config.get("domain", "ny"))[0]
nz = ast.literal_eval(config.get("domain", "nz"))[0]
z_origin = ast.literal_eval(config.get("domain", "z_origin"))[0]

y = np.arange(dy/2,dy*ny+dy/2,dy)
x = np.arange(dx/2,dx*nx+dx/2,dx)
z = np.arange(dz/2, dz*nz, dz) + z_origin
xu = x + np.gradient(x)/2
xu = xu[:-1]
yv = y + np.gradient(y)/2
yv = yv[:-1]
zw = z + np.gradient(z)/2
zw = zw[:-1]

## stretch factor for a vertically stretched grid
# set this to 1 if no streching required
dz_stretch_factor = ast.literal_eval(config.get("stretch", "dz_stretch_factor"))[0]

## Height level above which the grid is to be stretched vertically (in m)
dz_stretch_level = ast.literal_eval(config.get("stretch", "dz_stretch_level"))[0]

## allowed maximum vertical grid spacing (in m)
dz_max = ast.literal_eval(config.get("stretch", "dz_max"))[0]

if dz_stretch_factor>1.0:
    z, zw = calc_stretch(z, dz, zw, dz_stretch_level)

dz_soil = np.array(ast.literal_eval(config.get("soil", "dz_soil")))

wrf_path = ast.literal_eval(config.get("wrf", "wrf_path"))[0]
wrf_file = ast.literal_eval(config.get("wrf", "wrf_output"))

interp_mode = ast.literal_eval(config.get("wrf", "interp_mode"))[0]

start_year  = ast.literal_eval(config.get("wrf", "start_year"))[0]
start_month = ast.literal_eval(config.get("wrf", "start_month"))[0]
start_day   = ast.literal_eval(config.get("wrf", "start_day"))[0]
start_hour  = ast.literal_eval(config.get("wrf", "start_hour"))[0]

end_year  = ast.literal_eval(config.get("wrf", "end_year"))[0]
end_month = ast.literal_eval(config.get("wrf", "end_month"))[0]
end_day   = ast.literal_eval(config.get("wrf", "end_day"))[0]
end_hour  = ast.literal_eval(config.get("wrf", "end_hour"))[0]
dynamic_ts = ast.literal_eval(config.get("wrf", "dynamic_ts"))[0]


#-------------------------------------------------------------------------------
# Read WRF
#-------------------------------------------------------------------------------
## the input can be one wrf file, a list of files,
# or a string glob in the form "path/to/my/files/*.nc"
print("Reading WRF")
if len(wrf_file) == 1:
    wrf_files = sorted(glob(wrf_path+wrf_file[0]))
else:
    wrf_files = sorted([wrf_path+file for file in wrf_file ])

## use salem to read WRF
# remove duplicated timestamps
ds_wrf = xr.Dataset()
with salem.open_mf_wrf_dataset(wrf_files) as ds_raw:
    ## in case xtime is created as time dimension
    if len(ds_raw["time"])==1:
        ds_raw = ds_raw.isel(time=0)
        ds_raw = ds_raw.rename({"xtime": "time"})
    for variables in ds_raw.data_vars:
        ds_wrf[variables] = ds_raw[variables].drop_duplicates("time", keep="last")
    ds_wrf.attrs = ds_raw.attrs

del ds_raw

dt_start = datetime(start_year, start_month, start_day, start_hour,)
dt_end = datetime(end_year, end_month, end_day, end_hour,)

## check WRF temporal frequency; convert ns to s
wrf_ts = (ds_wrf["time"][1]-ds_wrf["time"][0]).data.astype("float64")* 1e-9

## temporal interpolation currently not supported in WRF4PALM
if dynamic_ts<wrf_ts:
    raise SystemExit(
    "Invalid timesteps given. Stopping..."
    )


## find how many timestamps to interpolate
num_ts = (dt_end - dt_start)/timedelta(seconds=dynamic_ts)
## generate a list of timestamps
all_ts = [dt_start+i*timedelta(seconds=dynamic_ts) for i in range(0,floor(num_ts)+1)]
## round up the end time index so that PALM doesn't crash
# when data of the final timestamp is not given
if floor(num_ts) != ceil(num_ts):
    all_ts.append(dt_end)

all_ts = np.array(all_ts).astype("datetime64[ns]")
## select required timestamps
ds_wrf = ds_wrf.sel(time=all_ts)
ds_wrf["PT"] = ds_wrf["T"] + 300
ds_wrf["PT"].attrs = ds_wrf["T"].attrs
ds_wrf["PT"].attrs["description"] = "potential temperature"
#-------------------------------------------------------------------------------
# Locate PALM domain in WRF
#-------------------------------------------------------------------------------
## find WRF map projection
map_proj = ds_wrf.MAP_PROJ

wrf_map_dict = {
                1: "lcc",
                2: "stere",
                3: "merc",
                6: "latlong",
}

if map_proj not in wrf_map_dict:
    raise SystemExit(
    "Incompatible WRF map projection, stopping..."
    )

wgs_proj = Proj(proj='latlong', datum='WGS84', ellips='sphere')

dx_wrf, dy_wrf = ds_wrf.DX, ds_wrf.DY
if map_proj == 6:
    wrf_proj = wgs_proj
    xx_wrf = ds_wrf.lon.data
    yy_wrf = ds_wrf.lat.data
else:
    wrf_proj = Proj(proj=wrf_map_dict[map_proj], # projection type
                    lat_1=ds_wrf.TRUELAT1, lat_2=ds_wrf.TRUELAT2,
                    lat_0=ds_wrf.MOAD_CEN_LAT, lon_0=ds_wrf.STAND_LON,
                    a=6370000, b=6370000) # The Earth is a perfect sphere in WRF

    # Easting and Northings of the domains center point
    trans_wgs2wrf = Transformer.from_proj(wgs_proj, wrf_proj)
    e, n = trans_wgs2wrf.transform(ds_wrf.CEN_LON, ds_wrf.CEN_LAT)
    # WRF Grid parameters
    nx_wrf, ny_wrf = ds_wrf.dims['west_east'], ds_wrf.dims['south_north']
    # Down left corner of the domain
    x0_wrf = -(nx_wrf-1) / 2. * dx_wrf + e
    y0_wrf = -(ny_wrf-1) / 2. * dy_wrf + n
    # 2d grid
    xx_wrf, yy_wrf = np.meshgrid(np.arange(nx_wrf) * dx_wrf + x0_wrf,
                                 np.arange(ny_wrf) * dy_wrf + y0_wrf)

## if no PALM projection is given by user,
#  then use WGS84 lat/lon and WRF projection to locate domain
# otherwise use the user specified projection
if len(palm_proj_code) == 0:
    palm_proj = wrf_proj
else:
    palm_proj = Proj(init = palm_proj_code)

trans_wrf2palm = Transformer.from_proj(wrf_proj, palm_proj)
lons_wrf,lats_wrf = trans_wrf2palm.transform(xx_wrf, yy_wrf)

west, east, south, north, centx, centy = domain_location(palm_proj, wgs_proj, centlat, centlon,
                                           dx, dy, nx, ny)

## write a cfg file for future reference
generate_cfg(case_name, dx, dy, dz, nx, ny, nz,
             west, east, south, north, centlat, centlon,z_origin)

# find indices of closest values
west_idx,east_idx,south_idx,north_idx = framing_2d_cartesian(lons_wrf,lats_wrf, west,east,south,north,dx_wrf,dy_wrf)
# in case negative longitudes are used
# these two lines may be redundant need further tests 27 Oct 2021
if east_idx-west_idx<0:
    east_idx, west_idx = west_idx, east_idx

# If PALM domain smaller than one WRF grid spacing
if north_idx-south_idx<1 or east_idx-west_idx<1:
    print(north_idx, south_idx,  east_idx, west_idx)
    raise SystemExit(
    "PALM domain size is smaller than one WRF grid cell size.\n"+
    "Please consider re-configure your PALM domain.\n"+
    "Stopping...\n"
    )
## drop data outside of PALM domain area
mask_sn = (ds_wrf.south_north>=ds_wrf.south_north[south_idx]) & (ds_wrf.south_north<=ds_wrf.south_north[north_idx])
mask_we = (ds_wrf.west_east>=ds_wrf.west_east[west_idx]) & (ds_wrf.west_east<=ds_wrf.west_east[east_idx])

ds_drop = ds_wrf.where(mask_sn & mask_we, drop=True)

#-------------------------------------------------------------------------------
# Read dynamic driver
#-------------------------------------------------------------------------------
dynstr = f'dynamic_files/{case_name}_dynamic_{start_year}_{start_month}_{start_day}_{start_hour}'
ds_dynamic = xr.open_dataset(dynstr)
#-------------------------------------------------------------------------------
# Plotting functions
#-------------------------------------------------------------------------------
def zcross(ds_drop, ds_dynamic, var, ts=0):
    # plot vertical cross sections
    if type(ts) is str:
        ts = np.argwhere(pd.to_datetime(all_ts)==pd.to_datetime(ts))[0,0]

    plt.figure(figsize = (10, 20))
    palm_var = var.lower()
    wrf_var = var.upper()
    if wrf_var=="QV":
        wrf_var = "QVAPOR"

    palm_bc_list = ["ls_forcing_left_", "ls_forcing_right_", "ls_forcing_south_", "ls_forcing_north_"]
    palm_var_list = [ls_str+palm_var for ls_str in palm_bc_list]

    wrf_ds_list = [ds_drop[wrf_var].isel(time=ts, west_east=0),
                   ds_drop[wrf_var].isel(time=ts, west_east=-1),
                   ds_drop[wrf_var].isel(time=ts, south_north=0),
                   ds_drop[wrf_var].isel(time=ts, south_north=-1)]

    vmin = ds_dynamic[palm_var_list[0]].isel(time=ts).min()
    vmax = ds_dynamic[palm_var_list[0]].isel(time=ts).max()

    plt.suptitle(all_ts[ts])
    i = 1
    for varplot in palm_var_list:
        plt.subplot(4,2, i)
        palm_ds = ds_dynamic[varplot].isel(time=ts)
        palm_ds.plot(vmin = vmin, vmax = vmax)
        plt.title(palm_var+" in PALM")

        i = i+2
    j = 2
    for ds in wrf_ds_list:
        plt.subplot(4,2, j)
        if j < 5:
            plt.pcolormesh(ds.south_north-ds.south_north[0],
                           ds_drop.Z.mean(dim=["time", "south_north", "west_east"]), ds,
                           vmin=vmin, vmax=vmax, shading ="auto")
            plt.xlabel("south-north (m)")
        else:
            plt.pcolormesh(ds.west_east-ds.west_east[0],
                           ds_drop.Z.mean(dim=["time", "south_north", "west_east"]), ds,
                           vmin=vmin, vmax=vmax, shading ="auto")
            plt.xlabel("west-east (m)")
        plt.colorbar()
        plt.title(wrf_var+" in WRF")
        plt.ylim([ds_dynamic.z[0], ds_dynamic.z[-1]])
        plt.ylabel("height (m)")

        j = j+2

    plt.tight_layout()
    plt.show()
def pr(ds_drop, ds_dynamic, var, ts=0):
    # plot vertical profiles
    if type(ts) is str:
        ts = np.argwhere(pd.to_datetime(all_ts)==pd.to_datetime(ts))[0,0]

    plt.figure(figsize = (16, 9))
    palm_var = var.lower()
    wrf_var = var.upper()
    if wrf_var=="QV":
        wrf_var = "QVAPOR"

    palm_bc_list = ["ls_forcing_left_", "ls_forcing_right_", "ls_forcing_south_", "ls_forcing_north_"]
    palm_var_list = [ls_str+palm_var for ls_str in palm_bc_list]

    wrf_ds_list = [ds_drop[wrf_var].isel(time=ts, west_east=0),
                   ds_drop[wrf_var].isel(time=ts, west_east=-1),
                   ds_drop[wrf_var].isel(time=ts, south_north=0),
                   ds_drop[wrf_var].isel(time=ts, south_north=-1)]

    vmin = ds_dynamic[palm_var_list[0]].isel(time=ts).min()
    vmax = ds_dynamic[palm_var_list[0]].isel(time=ts).max()

    plt.suptitle(all_ts[ts])
    for i, varplot in enumerate(palm_var_list):
        plt.subplot(1,4, i+1)
        palm_ds = ds_dynamic[varplot].isel(time=ts).mean(axis=1)
        ds = wrf_ds_list[i]
        plt.plot(ds.mean(axis=1), ds_drop.Z.mean(dim=["time", "south_north", "west_east"]),
                 'x-', label="WRF")
        if palm_var == "w":
            plt.plot(palm_ds, ds_dynamic.zw,"--", label ="PALM")
            plt.ylim([ds_dynamic.zw[0], ds_dynamic.zw[-1]])
        else:
            plt.plot(palm_ds, ds_dynamic.z,"--", label ="PALM")
            plt.ylim([ds_dynamic.z[0], ds_dynamic.z[-1]])
        plt.ylabel("z (m)")
        plt.xlabel(palm_var)
        plt.xlim([vmin*0.9, vmax*1.2])
        plt.title(varplot)
        plt.legend()
        plt.grid(alpha=0.7)
    plt.tight_layout()
    plt.show()

def ts(ds_drop, ds_dynamic, var, level=0):
    # plot time series
    plt.figure(figsize = (10, 16))
    palm_var = var.lower()
    wrf_var = var.upper()
    if wrf_var=="QV":
        wrf_var = "QVAPOR"

    palm_bc_list = ["ls_forcing_left_", "ls_forcing_right_",
                    "ls_forcing_south_", "ls_forcing_north_", "ls_forcing_top_"]
    palm_var_list = [ls_str+palm_var for ls_str in palm_bc_list]

    wrf_z = ds_drop.Z.mean(dim=["time", "south_north", "west_east"]).load().data
    _, wrf_top_idx = nearest_1d(wrf_z, ds_dynamic.z.data[-1])
    _, wrf_idx = nearest_1d(wrf_z, level)
    wrf_ds_list = [ds_drop[wrf_var].isel(bottom_top = wrf_idx, west_east=0),
                   ds_drop[wrf_var].isel(bottom_top = wrf_idx, west_east=-1),
                   ds_drop[wrf_var].isel(bottom_top = wrf_idx, south_north=0),
                   ds_drop[wrf_var].isel(bottom_top = wrf_idx, south_north=-1),
                   ds_drop[wrf_var].isel( bottom_top = wrf_top_idx)]


    for i, varplot in enumerate(palm_var_list):
        wrf_ds = wrf_ds_list[i]
        # i<4 -> lateral boundaries
        if i <4:
            if palm_var == "w":
                _, palm_idx = nearest_1d(ds_dynamic.zw.data, level)
                palm_ds = ds_dynamic[varplot].isel(zw=palm_idx).mean(axis=1)
            else:
                _, palm_idx = nearest_1d(ds_dynamic.z.data, level)
                palm_ds = ds_dynamic[varplot].isel(z=palm_idx).mean(axis=1)
            wrf_ds_plot = wrf_ds.mean(axis=1)
            plt_title = f"{varplot} at {level} m (nearest)"
        # top boundary
        else:
            palm_ds = ds_dynamic[varplot].mean(axis=(1,2))
            wrf_ds_plot = wrf_ds.mean(axis=(1,2))
            plt_title = f"{varplot}"
        plt.subplot(5,1, i+1)
        plt.plot(ds_drop.time, wrf_ds_plot, 'x-', label="WRF")
        plt.plot(ds_drop.time, palm_ds, "--", label="PALM")
        plt.title(plt_title)
        plt.xlabel("time")
        plt.legend()
        plt.grid(alpha=0.7)
    plt.tight_layout()
    plt.show()

plot_type = sys.argv[2]
plot_var = sys.argv[3]

if plot_type == "zcross":
    ts = input("Please enter the timestamp (yyyy-mm-dd-hh): ")
    zcross(ds_drop, ds_dynamic, plot_var, ts)
elif plot_type == "pr":
    ts = input("Please enter the timestamp (yyyy-mm-dd-hh): ")
    pr(ds_drop, ds_dynamic, plot_var, ts)
elif plot_type == "ts":
    level = input("Please enter the vertical level in m: ")
    ts(ds_drop, ds_dynamic, plot_var, float(level))
