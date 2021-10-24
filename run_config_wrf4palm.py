#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------
# WRF4PALM 
#--------------------------------------------------------------------------------
# Process data from WRF to PALM v6.0     
# Output of this script is the NetCDF dynamic driver for PALM following          
# PALM Input Data Standard (PIDS) v1.9
# [Update v1.1]- Oct 2021
# - use salem, xarray and multiprocessing
# - update namelist variables
# - add WRF projection configuration
# - modify geostrophic wind calculation
# 
# [Update v1.0.1] - 7 Jan 2021
# User input should be provided in namelist.dynamic
#                                                                              
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)                          
# Acknowledgement: The author would like to acknowledge Ricardo Faria for his     
# initial contribution of WRF2PALM https://github.com/ricardo88faria/WRF2PALM.   
#--------------------------------------------------------------------------------
import sys 
import time
import salem
import xarray as xr
from functools import partial
from pyproj import Proj, Transformer
import configparser
import ast
import numpy as np
from math import ceil, floor
from datetime import datetime, timedelta
from tqdm import tqdm
from functools import partial
from multiprocess import Pool
from dynamic_util.nearest import nearest_2d
from dynamic_util.loc_dom import calc_stretch, domain_location, generate_cfg
from dynamic_util.process_wrf import zinterp, multi_zinterp, process_top
from dynamic_util.geostrophic import calc_geostrophic_wind
from dynamic_util.surface_nan_solver import *
import warnings
## supress warnings
## switch to other actions if needed
warnings.filterwarnings("ignore", '.*pyproj.*')

start = datetime.now()
#--------------------------------------------------------------------------------
# Read user input namelist
#--------------------------------------------------------------------------------
settings_cfg = configparser.ConfigParser(inline_comment_prefixes='#')
config = configparser.RawConfigParser()
config.read(sys.argv[1])#"namelist.test")     
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
    z, zw = calc_stretch(z, dz)

dz_soil = np.array(ast.literal_eval(config.get("soil", "dz_soil")))
msoil_val = np.array(ast.literal_eval(config.get("soil", "msoil")))[0]


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
    for variables in ds_raw.data_vars:
        ds_wrf[variables] = ds_raw[variables].drop_duplicates("time", keep="last")
    ds_wrf.attrs = ds_raw.attrs

del ds_raw


#-------------------------------------------------------------------------------
# Find timestamps
#-------------------------------------------------------------------------------
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
# calculate timestamp in seconds
time_step_sec = ((dt_end-dt_start)).total_seconds()
times_sec = np.zeros(len(all_ts))
for t in range(0,len(all_ts)):
    times_sec[t] = (all_ts[t]-all_ts[0]).astype('float')*1e-9
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

if map_proj == 6:
    wrf_proj = wgs_proj
else:
    wrf_proj = Proj(proj=wrf_map_dict[map_proj], # projection type
                    lat_1=ds_wrf.TRUELAT1, lat_2=ds_wrf.TRUELAT2, 
                    lat_0=ds_wrf.MOAD_CEN_LAT, lon_0=ds_wrf.STAND_LON, 
                    a=6370000, b=6370000) # The Earth is a perfect sphere in WRF

# Easting and Northings of the domains center point
trans_wgs2wrf = Transformer.from_proj(wgs_proj, wrf_proj)
e, n = trans_wgs2wrf.transform(ds_wrf.CEN_LON, ds_wrf.CEN_LAT)
# e, n = transform(wgs_proj, wrf_proj, ds_wrf.CEN_LON, ds_wrf.CEN_LAT)
# WRF Grid parameters
dx_wrf, dy_wrf = ds_wrf.DX, ds_wrf.DY
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

west, east, south, north = domain_location(palm_proj, wgs_proj, centlat, centlon, dx, dy, nx, ny)

## write a cfg file for future reference
generate_cfg(case_name, dx, dy, dz, nx, ny, nz, west, east, south, north, centlat, centlon,z_origin)

# find indices of closest values
south_idx, north_idx = nearest_2d(lats_wrf, south)[1][0], nearest_2d(lats_wrf, north)[1][0]
west_idx, east_idx = nearest_2d(lons_wrf, west)[1][1], nearest_2d(lons_wrf, east)[1][1]  

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
ds_drop["pt"] = ds_drop["T"] + 300
ds_drop["pt"].attrs = ds_drop["T"].attrs
ds_drop["gph"] = (ds_drop["PH"] + ds_drop["PHB"])/9.81
ds_drop["gph"].attrs = ds_drop["PH"].attrs


#-------------------------------------------------------------------------------
# Horizontal interpolation
#-------------------------------------------------------------------------------
print("Start horizontal interpolation")
# assign new coordinates based on PALM
south_north_palm = ds_drop.south_north[0].data+y
west_east_palm = ds_drop.west_east[0].data+x
# staggered coordinates
south_north_v_palm = ds_drop.south_north[0].data+yv
west_east_u_palm = ds_drop.west_east[0].data+xu

# interpolation
ds_drop = ds_drop.assign_coords({"west_east_palm": west_east_palm,
                                 "south_north_palm": south_north_palm,
                                 "west_east_u_palm": west_east_u_palm,
                                 "south_north_v_palm": south_north_v_palm})
ds_interp = ds_drop.interp({"west_east": ds_drop.west_east_palm,}, method = interp_mode
                          ).interp({"south_north": ds_drop.south_north_palm}, method = interp_mode)
ds_interp_u = ds_drop.interp({"west_east": ds_drop.west_east_u_palm,}, method = interp_mode
                          ).interp({"south_north": ds_drop.south_north_palm}, method = interp_mode)
ds_interp_v = ds_drop.interp({"west_east": ds_drop.west_east_palm,}, method = interp_mode
                          ).interp({"south_north": ds_drop.south_north_v_palm}, method = interp_mode)

ds_interp = ds_interp.drop(["west_east", "south_north"]
                          ).rename({"west_east_palm": "west_east",
                                    "south_north_palm": "south_north"})

ds_interp_u = ds_interp_u.drop(["west_east", "south_north"]
                          ).rename({"west_east_u_palm": "west_east",
                                    "south_north_palm": "south_north"})

ds_interp_v = ds_interp_v.drop(["west_east", "south_north"]
                          ).rename({"west_east_palm": "west_east",
                                    "south_north_v_palm": "south_north"})

## get surface and soil fields
zs_wrf = ds_interp.ZS[0,:,0,0].load()
t2_wrf = ds_interp.T2.load()
u10_wrf = ds_interp_u.U10.load()
v10_wrf = ds_interp_v.V10.load()
qv2_wrf = ds_interp.Q2.load()
psfc_wrf = ds_interp.PSFC.load()
pt2_wrf = t2_wrf*((1000)/(psfc_wrf*0.01))**0.286

surface_var_dict = {"U": u10_wrf,
                   "V": v10_wrf,
                   "pt": pt2_wrf,
                   "QVAPOR": qv2_wrf,
                   "W": None}

#-------------------------------------------------------------------------------
# soil moisture and temperature
#------------------------------------------------------------------------------- 
print("Calculating soil temperature and moisture from WRF")

watermask = ds_interp["LANDMASK"].sel(time=dt_start).load().data == 0
landmask = ds_interp["LANDMASK"].sel(time=dt_start).load().data == 1
median_smois = [np.nanmedian(ds_interp["SMOIS"][0,izs,:,:].load().data[landmask]) for izs in range(0,len(zs_wrf))]
ds_interp["soil_layers"] = zs_wrf.load().data
tslb_wrf = ds_interp["TSLB"].sel(time=dt_start).load()
smois_wrf = ds_interp["SMOIS"].sel(time=dt_start).load()
deep_soil_wrf = ds_interp["TMN"].sel(time=dt_start)
deep_tsoil = deep_soil_wrf.where(landmask).mean().load().data
for izs in range(0,len(zs_wrf)):
    smois_wrf.isel(soil_layers=izs).data[watermask] = median_smois[izs]
    if smois_wrf.isel(soil_layers=izs).mean()== 0.0:
        smois_wrf.isel(soil_layers=izs).data[:,:] = msoil_val
    
init_tsoil = np.zeros((len(dz_soil), len(y), len(x)))
init_msoil = np.zeros((len(dz_soil), len(y), len(x)))
for iy in tqdm(range(0,len(y)),position=0, leave=True):
    for ix in range(0, len(x)):
        init_tsoil[:,iy,ix] = np.interp(dz_soil, zs_wrf.data, tslb_wrf[:,iy,ix])
        init_msoil[:,iy,ix] = np.interp(dz_soil, zs_wrf.data, smois_wrf[:,iy,ix])

#-------------------------------------------------------------------------------
# Vertical interpolation
#-------------------------------------------------------------------------------
print("Start vertical interpolation")
# create an empty dataset to store interpolated data
ds_interp["Z"] = ds_interp["Z"].load()
ds_we = ds_interp.isel(west_east=[0,-1])
ds_sn = ds_interp.isel(south_north=[0,-1])

ds_we_ustag = ds_interp_u.isel(west_east=[0,-1])
ds_we_vstag = ds_interp_v.isel(west_east=[0,-1])

ds_sn_ustag = ds_interp_u.isel(south_north=[0,-1])
ds_sn_vstag = ds_interp_v.isel(south_north=[0,-1])

varbc_list = ["W", "QVAPOR","pt","Z"]
for var in ds_we.data_vars:
    if var not in varbc_list:
        ds_we = ds_we.drop(var)
        ds_sn = ds_sn.drop(var)
    if var not in ["U", "Z"]:
        ds_we_ustag = ds_we_ustag.drop(var)
        ds_sn_ustag = ds_sn_ustag.drop(var)
    if var not in ["V", "Z"]:
        ds_we_vstag = ds_we_vstag.drop(var)
        ds_sn_vstag = ds_sn_vstag.drop(var)
        
ds_we = ds_we.load()
ds_sn = ds_sn.load()

ds_we_ustag = ds_we_ustag.load()
ds_sn_ustag = ds_sn_ustag.load()

ds_we_vstag = ds_we_vstag.load()
ds_sn_vstag = ds_sn_vstag.load()

ds_palm_we = xr.Dataset()
ds_palm_we = ds_palm_we.assign_coords({"x": x[:2],"y": y, "time":ds_interp.time.data, 
                                       "z": z, "yv": yv, "xu": xu[:2], "zw":zw})

ds_palm_sn = xr.Dataset()
ds_palm_sn = ds_palm_sn.assign_coords({"x": x,"y": y[:2], "time":ds_interp.time.data, 
                                       "z": z, "yv": yv[:2], "xu": xu, "zw":zw})

zeros_we = np.zeros((len(all_ts), len(z), len(y), len(x[:2]))) 
zeros_sn = np.zeros((len(all_ts), len(z), len(y[:2]), len(x))) 

# interpolation scalars
for varbc in ["QVAPOR","pt"]:
    ds_palm_we[varbc] = xr.DataArray(np.copy(zeros_we), dims=['time','z','y', 'x'])
    ds_palm_sn[varbc] = xr.DataArray(np.copy(zeros_sn), dims=['time','z','y', 'x'])
    print(f"Processing {varbc} for west and east boundaries")
    ds_palm_we[varbc] = multi_zinterp(max_pool, ds_we, varbc, z, ds_palm_we)
    print(f"Processing {varbc} for south and north boundaries")
    ds_palm_sn[varbc] = multi_zinterp(max_pool, ds_sn, varbc, z, ds_palm_sn)
    
# interpolate w
zeros_we_w = np.zeros((len(all_ts), len(zw), len(y), len(x[:2]))) 
zeros_sn_w = np.zeros((len(all_ts), len(zw), len(y[:2]), len(x))) 
ds_palm_we["W"] = xr.DataArray(np.copy(zeros_we_w), dims=['time','zw','y', 'x'])
ds_palm_sn["W"] = xr.DataArray(np.copy(zeros_sn_w), dims=['time','zw','y', 'x'])

print(f"Processing W for west and east boundaries")
ds_palm_we["W"] = multi_zinterp(max_pool, ds_we, "W", zw, ds_palm_we)
print(f"Processing W for south and north boundaries")
ds_palm_sn["W"] = multi_zinterp(max_pool, ds_sn, "W", zw, ds_palm_sn)

# interpolate u and v
zeros_we_u = np.zeros((len(all_ts), len(z), len(y), len(xu[:2]))) 
zeros_sn_u = np.zeros((len(all_ts), len(z), len(y[:2]), len(xu))) 
ds_palm_we["U"] = xr.DataArray(np.copy(zeros_we_u), dims=['time','z','y', 'xu'])
print(f"Processing U for west and east boundaries")
ds_palm_we["U"] = multi_zinterp(max_pool, ds_we_ustag, "U", z, ds_palm_we)

ds_palm_sn["U"] = xr.DataArray(np.copy(zeros_sn_u), dims=['time','z','y', 'xu'])
print(f"Processing U for south and north boundaries")
ds_palm_sn["U"] = multi_zinterp(max_pool, ds_sn_ustag, "U", z, ds_palm_sn)

zeros_we_v = np.zeros((len(all_ts), len(z), len(yv), len(x[:2]))) 
zeros_sn_v = np.zeros((len(all_ts), len(z), len(yv[:2]), len(x))) 
ds_palm_we["V"] = xr.DataArray(np.copy(zeros_we_v), dims=['time','z','yv', 'x'])
print(f"Processing V for west and east boundaries")
ds_palm_we["V"] = multi_zinterp(max_pool, ds_we_vstag, "V", z, ds_palm_we)

ds_palm_sn["V"] = xr.DataArray(np.copy(zeros_sn_v), dims=['time','z','yv', 'x'])
print(f"Processing V for south and north boundaries")
ds_palm_sn["V"] = multi_zinterp(max_pool, ds_sn_vstag, "V", z, ds_palm_sn)
#-------------------------------------------------------------------------------
# top boundary
#-------------------------------------------------------------------------------
print("Processing top boundary conditions...")
u_top = np.zeros((len(all_ts), len(y), len(xu)))
v_top = np.zeros((len(all_ts), len(yv), len(x)))
w_top = np.zeros((len(all_ts), len(y), len(x)))
qv_top = np.zeros((len(all_ts), len(y), len(x)))
pt_top = np.zeros((len(all_ts), len(y), len(x)))

for var in ds_interp.data_vars:
    if var not in varbc_list:
        ds_interp = ds_interp.drop(var)
    if var not in ["U", "Z"]:
        ds_interp_u = ds_interp_u.drop(var)
    if var not in ["V", "Z"]:
        ds_interp_v = ds_interp_v.drop(var)

ds_interp = ds_interp.load()
ds_interp_u = ds_interp_u.load()
ds_interp_v = ds_interp_v.load()

        
top_dict = {"U": (ds_interp_u, u_top, z),
            "V": (ds_interp_v, v_top, z),
            "pt": (ds_interp, pt_top, z),
            "QVAPOR": (ds_interp, qv_top, z),
            "W": (ds_interp, w_top, zw)}

with Pool(max_pool) as p:
        pool_outputs = list(
            tqdm(
                p.imap(partial(process_top, all_ts,top_dict), top_dict.keys()), total=len(top_dict.keys()),
                position=0, leave=True
            )
        )
p.join()
    ## convert dictionary back to dataset
pool_dict = dict(pool_outputs)
u_top = pool_dict["U"]
v_top = pool_dict["V"]
w_top = pool_dict["W"]
qv_top = pool_dict["QVAPOR"]
pt_top = pool_dict["pt"]
#-------------------------------------------------------------------------------
# Geostrophic wind estimation
#-------------------------------------------------------------------------------
print("Geostrophic wind estimation...")
lat_geostr = ds_drop.lat[:,0]
dx_wrf = ds_drop.DX
dy_wrf = ds_drop.DY
gph = ds_drop.gph
gph = gph.load()
ds_geostr = xr.Dataset()
ds_geostr = ds_geostr.assign_coords({"time":ds_drop.time.data, 
                                     "z": ds_drop["Z"].mean(("time", "south_north", "west_east")).data})
ds_geostr["ug"] = xr.DataArray(np.zeros((len(all_ts),len(gph.bottom_top.data))),
                               dims=['time','z'])
ds_geostr["vg"] = xr.DataArray(np.zeros((len(all_ts),len(gph.bottom_top.data))),
                               dims=['time','z'])

for ts in tqdm(range(0,len(all_ts)), total=len(all_ts), position=0, leave=True):
    for levels in gph.bottom_top.data:
        ds_geostr["ug"][ts,levels], ds_geostr["vg"][ts,levels] = calc_geostrophic_wind(
        gph[ts,levels, :,:].data, lat_geostr.data, dy_wrf, dx_wrf)


# interpolate to PALM vertical levels
ds_geostr = ds_geostr.interp({"z": z}) 
        
        
#-------------------------------------------------------------------------------
# surface NaNs
#------------------------------------------------------------------------------- 
print("Resolving surface NaNs...")
# apply multiprocessing
with Pool(max_pool) as p:
    pool_outputs = list(
        tqdm(
            p.imap(partial(solve_surface,all_ts, ds_palm_we, ds_palm_sn, surface_var_dict),surface_var_dict.keys()), 
            total=len(surface_var_dict.keys()),position=0, leave=True
        )
    )
p.join()
pool_dict = dict(pool_outputs)
for var in surface_var_dict.keys():
    ds_palm_we[var]= pool_dict[var][0]
    ds_palm_sn[var]= pool_dict[var][1]
# near surface geostrophic wind
for t in range(0,len(all_ts)):
    ds_geostr["ug"][t,:] =  surface_nan_w(ds_geostr["ug"][t,:].data)
    ds_geostr["vg"][t,:] =  surface_nan_w(ds_geostr["vg"][t,:].data) 

#-------------------------------------------------------------------------------
# calculate initial profiles
#------------------------------------------------------------------------------- 
ds_drop["bottom_top"] = ds_drop["Z"].mean(("time", "south_north", "west_east")).data

u_init = ds_drop["U"].sel(time=dt_start).mean(
    dim=["south_north", "west_east"]).interp(
    {"bottom_top": z}, method = interp_mode)
v_init = ds_drop["V"].sel(time=dt_start).mean(
    dim=["south_north", "west_east"]).interp(
    {"bottom_top": z}, method = interp_mode)
# stagger w 
w_init = ds_drop["W"].sel(time=dt_start).mean(
    dim=["south_north", "west_east"]).interp(
    {"bottom_top": zw}, method = interp_mode)
qv_init = ds_drop["QVAPOR"].sel(time=dt_start).mean(
    dim=["south_north", "west_east"]).interp(
    {"bottom_top": z}, method = interp_mode)
pt_init = ds_drop["pt"].sel(time=dt_start).mean(
    dim=["south_north", "west_east"]).interp(
    {"bottom_top": z}, method = interp_mode)

u_init = surface_nan_uv(u_init.load().data, z, u10_wrf.sel(time=dt_start).mean(
                        dim=["south_north", "west_east"]).data)

v_init = surface_nan_uv(v_init.load().data, z, v10_wrf.sel(time=dt_start).mean(
                        dim=["south_north", "west_east"]).data)
w_init = surface_nan_w(w_init.load().data)
qv_init = surface_nan_s(qv_init.load().data, z, qv2_wrf.sel(time=dt_start).mean(
                        dim=["south_north", "west_east"]).data)
pt_init = surface_nan_s(pt_init.load().data, z, pt2_wrf.sel(time=dt_start).mean(
                        dim=["south_north", "west_east"]).data)

surface_pres = psfc_wrf[:, :,:].mean(dim=["south_north", "west_east"]).load()


#-------------------------------------------------------------------------------
# soil moisture and temperature
#------------------------------------------------------------------------------- 
nc_output_name = f'dynamic_files/{case_name}_dynamic_{start_year}_{start_month}_{start_day}_{start_hour}'
print('Writing NetCDF file',flush=True)
nc_output = xr.Dataset()
res_origin = str(dx) + 'x' + str(dy) + ' m'
nc_output.attrs['description'] = f'Contains dynamic data from WRF mesoscale. WRF output file: {wrf_file}'
nc_output.attrs['author'] = 'Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)'
nc_output.attrs['history'] = 'Created at ' + time.ctime(time.time())
nc_output.attrs['source']= 'netCDF4 python'
nc_output.attrs['origin_lat'] = float(centlat)
nc_output.attrs['origin_lon'] = float(centlon)
nc_output.attrs['z'] = float(0)
nc_output.attrs['x'] = float(0)
nc_output.attrs['y'] = float(0)
nc_output.attrs['rotation_angle'] = float(0)
nc_output.attrs['origin_time'] =  str(all_ts[0]) + ' UTC'
nc_output.attrs['end_time'] =  str(all_ts[-1]) + ' UTC'


nc_output['x'] = xr.DataArray(x, dims=['x'], attrs={'units':'m'})
nc_output['y'] = xr.DataArray(y, dims=['y'], attrs={'units':'m'})
nc_output['z'] = xr.DataArray(z-z_origin, dims=['z'], attrs={'units':'m'})
nc_output['zsoil'] = xr.DataArray(dz_soil, dims=['zsoil'], attrs={'units':'m'})
nc_output['xu'] = xr.DataArray(xu, dims=['xu'], attrs={'units':'m'})
nc_output['yv'] = xr.DataArray(yv, dims=['yv'], attrs={'units':'m'})
nc_output['zw'] = xr.DataArray(zw-z_origin, dims=['zw'], attrs={'units':'m'})
nc_output['time'] = xr.DataArray(times_sec, dims=['time'], attrs={'units':'seconds'})


nc_output.to_netcdf(nc_output_name)
nc_output['init_soil_m'] = xr.DataArray(init_msoil, dims=['zsoil','y','x'], 
         attrs={'units':'m^3/m^3','lod':np.int32(2), 'source':'WRF', 'long_name':'volumetric soil moisture (m^3/m^3)'}) 
nc_output['init_soil_t'] = xr.DataArray(init_tsoil, dims=['zsoil','y','x'], 
         attrs={'units':'K', 'lod':np.int32(2), 'source':'WRF', 'long_name':'soil temperature (K)'}) 

# output boundary conditions to PALM input
# directions: 0 west, 1 east
#             0 south, 1 north

nc_output['init_atmosphere_pt'] = xr.DataArray(pt_init,dims=['z'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_pt'] = xr.DataArray(ds_palm_we["pt"][:,:,:,0].data,dims=['time', 'z', 'y'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_pt'] = xr.DataArray(ds_palm_we["pt"][:,:,:,-1].data,dims=['time', 'z', 'y'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_pt'] = xr.DataArray(ds_palm_sn["pt"][:,:,0,:].data,dims=['time', 'z', 'x'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_pt'] = xr.DataArray(ds_palm_sn["pt"][:,:,-1,:].data,dims=['time', 'z', 'x'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
## top
nc_output['ls_forcing_top_pt'] = xr.DataArray(pt_top[:,:,:],dims=['time', 'y', 'x'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_qv'] = xr.DataArray(qv_init,dims=['z'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_qv'] = xr.DataArray(ds_palm_we["QVAPOR"][:,:,:,0].data,dims=['time', 'z', 'y'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_qv'] = xr.DataArray(ds_palm_we["QVAPOR"][:,:,:,-1].data,dims=['time', 'z', 'y'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_qv'] = xr.DataArray(ds_palm_sn["QVAPOR"][:,:,0,:].data,dims=['time', 'z', 'x'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_qv'] = xr.DataArray(ds_palm_sn["QVAPOR"][:,:,-1,:].data,dims=['time', 'z', 'x'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_qv'] = xr.DataArray(qv_top[:,:,:],dims=['time', 'y', 'x'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_u'] = xr.DataArray(u_init,dims=['z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_u'] = xr.DataArray(ds_palm_we["U"][:,:,:,0].data,dims=['time', 'z', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_u'] = xr.DataArray(ds_palm_we["U"][:,:,:,-1].data,dims=['time', 'z', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_u'] = xr.DataArray(ds_palm_sn["U"][:,:,0,:].data,dims=['time', 'z', 'xu'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_u'] = xr.DataArray(ds_palm_sn["U"][:,:,-1,:].data,dims=['time', 'z', 'xu'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_u'] = xr.DataArray(u_top[:,:,:],dims=['time', 'y', 'xu'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_v'] = xr.DataArray(v_init,dims=['z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_v'] = xr.DataArray(ds_palm_we["V"][:,:,:,0].data,dims=['time', 'z', 'yv'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_v'] = xr.DataArray(ds_palm_we["V"][:,:,:,-1].data,dims=['time', 'z', 'yv'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_v'] = xr.DataArray(ds_palm_sn["V"][:,:,0,:].data,dims=['time', 'z', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_v'] = xr.DataArray(ds_palm_sn["V"][:,:,-1,:].data,dims=['time', 'z', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_v'] = xr.DataArray(v_top[:,:,:],dims=['time', 'yv', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_w'] = xr.DataArray(w_init,dims=['zw'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_w'] = xr.DataArray(ds_palm_we["W"][:,:,:,0].data,dims=['time', 'zw', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_w'] = xr.DataArray(ds_palm_we["W"][:,:,:,-1].data,dims=['time', 'zw', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_w'] = xr.DataArray(ds_palm_sn["W"][:,:,0,:].data,dims=['time', 'zw', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_w'] = xr.DataArray(ds_palm_sn["W"][:,:,-1,:].data,dims=['time', 'zw', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_w'] = xr.DataArray(w_top[:,:,:],dims=['time', 'y', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})

nc_output['surface_forcing_surface_pressure'] = xr.DataArray(surface_pres.data, dims=['time'],
         attrs={'units':'Pa', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})


nc_output['ls_forcing_ug'] = xr.DataArray(ds_geostr["ug"].data,dims=['time','z'],
         attrs={'units':'m/s', 'long_name':'u wind component geostrophic', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_vg'] = xr.DataArray(ds_geostr["vg"].data,dims=['time','z'],
         attrs={'units':'m/s', 'long_name':'v wind component geostrophic', 'source':'WRF', 'res_origin':res_origin})



for var in nc_output.data_vars:
    encoding = {var: {'dtype': 'float32', '_FillValue': -9999, 'zlib':True}}
    nc_output[var].to_netcdf(nc_output_name, encoding=encoding, mode='a')


print('Add to your *_p3d file: ' + '\n soil_temperature = ' + 
              str([value for value in init_tsoil.mean(axis=(1,2))]) +
      '\n soil_moisture = ' + str([value for value in init_msoil.mean(axis=(1,2))]) 
        + '\n deep_soil_temperature = ' + str(deep_tsoil)+'\n')

with open('cfg_files/'+ case_name + '.cfg', "a") as cfg:
    cfg.write('Add to your *_p3d file: ' + '\n soil_temperature = ' + 
              str([value for value in init_tsoil.mean(axis=(1,2))]) +
      '\n soil_moisture = ' + str([value for value in init_msoil.mean(axis=(1,2))]) 
        + '\n deep_soil_temperature = ' + str(deep_tsoil)+'\n')




end = datetime.now()
print('PALM dynamic input file is ready. Script duration: {}'.format(end - start))
print('Start time: '+str(all_ts[0]))
print('End time: '+str(all_ts[-1]))
print('Time step: '+str(times_sec[1]-times_sec[0])+' seconds')
