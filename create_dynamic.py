#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
 Initially created on Mon Jul  2 12:08:53 2018 by @author: ricardofaria

 Script to create Large-scale forcing data (boundary conditions) netcdf input 
 for PALM following PALM Input Data Standard (PIDS) v1.9.
 – contains all dynamic data need to force PALM.
 Contains dynamic information for the run, such as time­dependent boundary 
 conditions and the initial state of the atmosphere.

#############################################################################
 Updated by Dongqi Lin
    Add the initial state of atmosphere parameters in netCDF dynamic driver
    Change interpolation method of boundary conditions
    Change geostropihc wind calculation
    Add start/end time index so that users can choose the time step to interpolate
    Time indexing can only apply to one WRF file
    No vertical streching is applied.
'''


###############################################################################
import numpy as np
from netCDF4 import Dataset, num2date
import wrf
from wrf import getvar, destagger, interplevel
import time

import pandas as pd
import xarray as xr
from tqdm import tqdm
from util import res_grid_change as rgc
from util import geostrophic
from util import nearest
from datetime import datetime
from util import surface_nan_solver as sns


start = datetime.now()

case_name = 'chch_50m_example'

wrf_file = 'wrfout_d03_2017-06-15_chch'

interp_mode = 'linear'

# !!! give start and end time to interpolate WRF output here !!!
# this depends on
# 1) WRF output time frequency
# 2) the desired PALM input update frequncy

dt_start = datetime(2017, 6, 16, 0,)
dt_end = datetime(2017, 6, 17, 0,)
interval = 2
ts = '2hour'

# layers for soil temperature and moisture calculation
# this shall be changed depending on different cases
dz_soil = np.array([0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86])

cfg = pd.read_csv('cfg_input/'+case_name + '.cfg')
dx = cfg.dx.values[0]
dy = cfg.dy.values[0]
dz = cfg.dz.values[0]
nx = cfg.nx.values[0]
ny = cfg.ny.values[0]
nz = cfg.nz.values[0]
north = cfg.north.values[0]
south = cfg.south.values[0]
east = cfg.east.values[0]
west = cfg.west.values[0]
lat_palm = cfg.midlat[0]
lon_palm = cfg.midlon[0]

y = np.arange(dy/2,dy*ny+1,dy)
x = np.arange(dx/2,dx*nx+1,dx)
z = np.arange(dz/2, dz*nz, dz)
xu = x + np.gradient(x)/2
xu = xu[:-1]
yv = y + np.gradient(y)/2
yv = yv[:-1]
zw = z + np.gradient(z)/2
zw = zw[:-1]

###############################################################################
##---------------------------- read WRF variables ---------------------------##
###############################################################################
times = []


print(f'Loading WRF netCDF: {wrf_file}' )
ds_wrf = xr.open_dataset(wrf_file)
nc_wrf = Dataset(wrf_file, 'r')

lat_s = ds_wrf['XLAT'][0,:,0].data
lon_s = ds_wrf['XLONG'][0,0,:].data
south_idx, north_idx = nearest.nearest(lat_s, south)[1], nearest.nearest(lat_s, north)[1]
west_idx, east_idx = nearest.nearest(lon_s, west)[1], nearest.nearest(lon_s, east)[1]

lat_v = ds_wrf['XLAT_V'][0,:,0].data
southv_idx, northv_idx = nearest.nearest(lat_v, south)[1], nearest.nearest(lat_v, north)[1]

lon_u = ds_wrf['XLONG_U'][0,0,:].data
westu_idx, eastu_idx = nearest.nearest(lon_u, west)[1], nearest.nearest(lon_u, east)[1]

lat_wrf = ds_wrf['XLAT'][0,south_idx:north_idx,0].data
lon_wrf = ds_wrf['XLONG'][0,0,west_idx:east_idx].data
# If PALM domain smaller than one WRF grid
if north_idx-south_idx<=1:
    north_idx = south_idx+2
if east_idx-west_idx<=1:
    east_idx = west_idx+2
if northv_idx-southv_idx<=1:
    northv_idx = southv_idx+2
if eastu_idx-westu_idx<=1:
    eastu_idx = westu_idx+2

# depth of soil layer
zs = ds_wrf['ZS'][0, :].data
# thickness of soil layer
dzs = ds_wrf['DZS'][0, :].data
# landmask - 1 is land, 0 is water
landmask = ds_wrf['LANDMASK'][0, south_idx:north_idx, west_idx:east_idx].data
# TMN - soil temperature at lower boundary
tmn = ds_wrf['TMN'][0, south_idx:north_idx, west_idx:east_idx].data

# Staggered ASL taking directly from WRF
PH = ds_wrf['PH']
PHB = ds_wrf['PHB']
H_stag = (PH + PHB) / 9.81
H = destagger(H_stag, stagger_dim=1)
wrf_time = nc_wrf.variables['XTIME']

time_var = num2date(nc_wrf.variables['XTIME'][:], nc_wrf.variables['XTIME'].units)
time_var = np.array(time_var).astype('datetime64[s]')
for i in range(0,time_var.shape[0]):
    if time_var[i] == dt_start:
        start_idx = i
        
    if time_var[i] == dt_end:
        end_idx = i

time_idx = np.arange(start_idx,end_idx+1,interval)

# round up the end time index so that PALM doesn't crash when the final time step is not given
input_lag = (dt_end-dt_start).total_seconds()
tmp_lag = (time_var[time_idx[-1]]-time_var[time_idx[0]]).astype('float')
if input_lag-tmp_lag > 0:
     time_idx= np.append(time_idx,end_idx)

## prepare arrays
pt = np.empty((time_idx.shape[0],H.shape[1],north_idx-south_idx,east_idx-west_idx))
pres = np.empty((time_idx.shape[0],H.shape[1],north_idx-south_idx,east_idx-west_idx))
tk = np.empty((time_idx.shape[0],H.shape[1],north_idx-south_idx,east_idx-west_idx))
tslb = np.empty((time_idx.shape[0],zs.shape[0],north_idx-south_idx,east_idx-west_idx))
smois = np.empty((time_idx.shape[0],zs.shape[0],north_idx-south_idx,east_idx-west_idx))
qv = np.empty((time_idx.shape[0],H.shape[1],north_idx-south_idx,east_idx-west_idx))
u = np.empty((time_idx.shape[0], H.shape[1], north_idx - south_idx, eastu_idx - westu_idx))
v = np.empty((time_idx.shape[0], H.shape[1], northv_idx - southv_idx, east_idx - west_idx))
w = np.empty((time_idx.shape[0], H_stag.shape[1], north_idx - south_idx, east_idx - west_idx))
zstag_wrf = np.empty_like(w)
z_wrf = np.empty_like(pt)
z_wrf_u = np.empty_like(u)
z_wrf_v = np.empty_like(v)



for t, wrf_t in enumerate(time_idx):
    pt[t,:,:,:] = getvar(nc_wrf, 'theta', timeidx = wrf_t, units='K')[:,south_idx:north_idx, west_idx:east_idx]
    pres[t,:,:,:] = getvar(nc_wrf, 'pres', timeidx = wrf_t, units='Pa')[:,south_idx:north_idx, west_idx:east_idx]
    tk[t,:,:,:] = getvar(nc_wrf, 'tk', timeidx = wrf_t)[:,south_idx:north_idx, west_idx:east_idx]
    # soil tempearture
    tslb[t,:,:,:] = nc_wrf.variables['TSLB'][wrf_t, :, south_idx:north_idx, west_idx:east_idx].data
    # soil moisture
    smois[t,:,:,:] = nc_wrf.variables['SMOIS'][wrf_t, :, south_idx:north_idx, west_idx:east_idx].data
    qv[t,:,:,:] = nc_wrf.variables['QVAPOR'][wrf_t, :, south_idx:north_idx, west_idx:east_idx].data
    # velocity field
    u[t,:,:,:] = ds_wrf['U'][wrf_t, :, south_idx:north_idx, westu_idx:eastu_idx].data
    v[t,:,:,:] = ds_wrf['V'][wrf_t, :, southv_idx:northv_idx, west_idx:east_idx].data
    w[t,:,:,:] = ds_wrf['W'][wrf_t, :, south_idx:north_idx, west_idx:east_idx].data
    zstag_wrf[t,:,:,:] = H_stag[wrf_t,:, south_idx:north_idx, west_idx:east_idx]
    z_wrf[t,:,:,:] = H[wrf_t,:, south_idx:north_idx, west_idx:east_idx]
    z_wrf_u[t,:,:,:] = H[wrf_t,:, south_idx:north_idx, westu_idx:eastu_idx]
    z_wrf_v[t,:,:,:] = H[wrf_t,:, southv_idx:northv_idx, west_idx:east_idx]
    times = np.append(times, num2date(wrf_time[wrf_t],wrf_time.units))

nc_wrf.close()

print("WRF output reading done.",flush=True)


time_step_sec = ((times[1]-times[0])).total_seconds()
times_sec = np.zeros(time_idx.shape[0])
for t, wrf_t in enumerate(time_idx):
    times_sec[t] = (time_var[wrf_t]-time_var[time_idx[0]]).astype('float')



def search_nan(var,t,var_type):
    if np.argwhere(np.isnan(var[t,0,:,:])).size > 0:
        for y_idx in range(0,var.shape[2]):
            for x_idx in range(0,var.shape[3]):
                if var_type == 'uv':
                    var[t,:,y_idx,x_idx] = sns.surface_nan_uv(var[t,:,y_idx,x_idx],z)
                elif var_type == 's':
                    var[t,:,y_idx,x_idx] = sns.surface_nan_s(var[t,:,y_idx,x_idx])
                elif var_type == 'w':
                    var[t,:,y_idx,x_idx] = sns.surface_nan_w(var[t,:,y_idx,x_idx])
                else:
                    print("wrong type given")
    return(var[t,:,:,:])

  

###############################################################################
##--------------------------- horizontal interpolation ----------------------##
###############################################################################
print("Interpolating horizontal fields...",flush=True)

# Arrays for horizontal interpolation
u_int = np.empty((times.shape[0], u.shape[1],y.shape[0],xu.shape[0]))
v_int = np.empty((times.shape[0], v.shape[1],yv.shape[0],x.shape[0]))
w_int = np.empty((times.shape[0], w.shape[1],y.shape[0],x.shape[0]))
qv_int = np.empty((times.shape[0], qv.shape[1],y.shape[0],x.shape[0]))
pt_int = np.empty((times.shape[0], pt.shape[1],y.shape[0],x.shape[0]))
pres_int = np.empty((times.shape[0], pres.shape[1],y.shape[0],x.shape[0]))
tk_int = np.empty((times.shape[0], tk.shape[1],y.shape[0],x.shape[0]))
z_wrf_int = np.empty((times.shape[0],u.shape[1],y.shape[0],x.shape[0]))
z_wrf_int_u = np.empty((times.shape[0],u.shape[1],y.shape[0],xu.shape[0]))
z_wrf_int_v = np.empty((times.shape[0],v.shape[1],yv.shape[0],x.shape[0]))
zstag_wrf_int = np.empty((times.shape[0],w.shape[1],y.shape[0],x.shape[0]))

for t in tqdm(range(u.shape[0]),ascii=True):
    for z_idx in range(0,u.shape[1]):
        u_int[t,z_idx,:,:] = rgc.grid_points_change(u[t,z_idx,:,:], xu.shape[0], y.shape[0], interp_mode)
        v_int[t,z_idx,:,:] = rgc.grid_points_change(v[t,z_idx,:,:], x.shape[0], yv.shape[0], interp_mode)
        qv_int[t,z_idx,:,:] = rgc.grid_points_change(qv[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        pt_int[t,z_idx,:,:] = rgc.grid_points_change(pt[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        pres_int[t,z_idx,:,:] = rgc.grid_points_change(pres[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        tk_int[t,z_idx,:,:] = rgc.grid_points_change(tk[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        z_wrf_int[t,z_idx,:,:] =  rgc.grid_points_change(z_wrf[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        z_wrf_int_u[t,z_idx,:,:] = rgc.grid_points_change(z_wrf_u[t,z_idx,:,:], xu.shape[0], y.shape[0], interp_mode)
        z_wrf_int_v[t,z_idx,:,:] = rgc.grid_points_change(z_wrf_v[t,z_idx,:,:], x.shape[0], yv.shape[0], interp_mode)
    for z_idx in range(0,w.shape[1]):
        w_int[t,z_idx,:,:] = rgc.grid_points_change(w[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
        zstag_wrf_int[t,z_idx,:,:] =  rgc.grid_points_change(zstag_wrf[t,z_idx,:,:], x.shape[0], y.shape[0], interp_mode)
print("Done.", flush=True)

###############################################################################
##--------------------------- vertical interpolation-------------------------##
###############################################################################
     
# Interpolation: unstaggered vetical levels
u_tmp = np.empty((times.shape[0], z.shape[0], u_int.shape[2], u_int.shape[3]))
v_tmp = np.empty((times.shape[0], z.shape[0],v_int.shape[2], v_int.shape[3]))
w_tmp = np.empty((times.shape[0], zw.shape[0],y.shape[0],x.shape[0]))
qv_tmp = np.empty((times.shape[0], z.shape[0],qv_int.shape[2], qv_int.shape[3]))
pt_tmp = np.empty((times.shape[0], z.shape[0],pt_int.shape[2], pt_int.shape[3]))
pres_tmp = np.empty((times.shape[0], z.shape[0],pres_int.shape[2], pres_int.shape[3]))
tk_tmp = np.empty((times.shape[0], z.shape[0],tk_int.shape[2], tk_int.shape[3]))

for l_idx, l in tqdm(enumerate(z), desc="Interpolating unstaggered vertical levels"):
    for t in range(0, times.shape[0]):
        qv_tmp[t, int(l_idx), :, :] = interplevel(qv_int[t, :, :, :], z_wrf_int[t,:,:,:], l).data
        pt_tmp[t, int(l_idx), :, :] = interplevel(pt_int[t, :, :, :], z_wrf_int[t,:,:,:], l).data
        u_tmp[t, int(l_idx), :, :] = interplevel(u_int[t, :, :, :], z_wrf_int_u[t,:,:,:], l).data
        v_tmp[t, int(l_idx), :, :] = interplevel(v_int[t, :, :, :], z_wrf_int_v[t,:,:,:], l).data
        pres_tmp[t, int(l_idx), :, :] = interplevel(pres_int[t, :, :, :], z_wrf_int[t,:,:,:], l).data
        tk_tmp[t, int(l_idx), :, :] = interplevel(tk_int[t, :, :, :], z_wrf_int[t,:,:,:], l).data

for lstag_idx, lstag in tqdm(enumerate(zw), desc="Interpolating staggered vertical levels"):    
    for t in range(0,times.shape[0]):
        w_tmp[t, int(lstag_idx),:, :] = interplevel(w_int[t,:,:,:], zstag_wrf_int[t,:,:,:], lstag).data

print(flush=True)
print('Vertical interpolation done.',flush=True)


# calculate geostrophic winds at every levels
# latitudes and longitudes are still required here
def rolling_mean(var, window):
    roll_mean = []
    for i in range(0, var.shape[0] - window, window):
        roll_mean.append(np.nansum(var[i:i + window]) / window)

    return (np.array(roll_mean))

lat_wrf_f = rgc.grid_points_change_1darray(lat_wrf,y.shape[0])
lon_wrf_f = rgc.grid_points_change_1darray(lon_wrf,x.shape[0])

geo_wind_u = np.zeros((pres_tmp.shape[0], pres_tmp.shape[1]))
geo_wind_v = np.zeros((pres_tmp.shape[0], pres_tmp.shape[1]))
geo_wind_u_f = np.zeros((u.shape[0], z.shape[0]))
geo_wind_v_f = np.zeros((u.shape[0], z.shape[0]))
for t in tqdm(range(pres_tmp.shape[0]),ascii=True, desc="Calculating geostropihc winds"):
    for h in range(0, pres_tmp.shape[1]):
        geo_wind = geostrophic.geostr(pres_tmp[t, h, :, :], tk_tmp[t, h, :, :], lat_wrf_f[:], lon_wrf_f[:], dy, dx)
        geo_wind_u[t, h] = geo_wind[0]
        geo_wind_v[t, h] = geo_wind[1]
        
    # "smooth" the geostrophic winds after calculation by taking rolling mean
    geo_wind_u_f[t, :] = rgc.grid_points_change_1darray(rolling_mean(geo_wind_u[t, :], 10), z.shape[0])
    geo_wind_v_f[t, :] = rgc.grid_points_change_1darray(rolling_mean(geo_wind_v[t, :], 10), z.shape[0])

print(flush=True)
print("Geostrophic wind calculation done.",flush=True)



for t in tqdm(range(0,u_tmp.shape[0]),ascii=True,desc = "Resolving surface NaNs"):
    u_tmp[t,:,:,:] = search_nan(u_tmp,t,'uv')
    v_tmp[t,:,:,:] = search_nan(v_tmp,t,'uv')
    pt_tmp[t,:,:,:] = search_nan(pt_tmp,t,'s')
    qv_tmp[t,:,:,:] = search_nan(qv_tmp,t,'s')
    pres_tmp[t,:,:,:] = search_nan(pres_tmp,t,'s')
    tk_tmp[t,:,:,:] = search_nan(tk_tmp,t,'s')
    w_tmp[t,:,:,:] = search_nan(w_tmp,t,'w')

print(flush='True')

# Genearte initial profiles
u_init = np.zeros(z.shape[0])
v_init = np.zeros(z.shape[0])
w_init = np.zeros(zw.shape[0])
qv_init = np.zeros(z.shape[0])
pt_init = np.zeros(z.shape[0])


for i in range(0,z.shape[0]):
    u_init[i] = np.nanmean(u_tmp[0,i,:,:].reshape(1,-1),axis=1)
    v_init[i] = np.nanmean(v_tmp[0,i,:,:].reshape(1,-1),axis=1)
    qv_init[i] = np.nanmean(qv_tmp[0,i,:,:].reshape(1,-1),axis=1)
    pt_init[i] = np.nanmean(pt_tmp[0,i,:,:].reshape(1,-1),axis=1)
for i in range(0, zw.shape[0]):
    w_init[i] = np.nanmean(w_tmp[0,i,:,:].reshape(1,-1),axis=1)


pres_surf = pres_tmp[:,1,:,:]
pres_init = []
for t in range(0, times.shape[0]):
    pres_init.append(np.nanmean(pres_surf[t,:,:]))

pres_init = np.array(pres_init)

# Soil temperature and moisture calculation
print('Calculating soil temperature and moisture from WRF', flush=True)
init_soil_t = np.zeros((dz_soil.shape[0], smois.shape[2], smois.shape[3]))
init_soil_m = np.zeros((dz_soil.shape[0], smois.shape[2], smois.shape[3]))
smois_land = np.zeros((zs.shape[0], smois.shape[2], smois.shape[3]))
init_soil_tmn = np.nanmean(np.ma.masked_where(landmask == 0, tmn))

def calc_soil_moisture(smois_lvl):
    # Function to exclude soli moisture==1 at water bodies
    # This is to avoid dismatch between WRF and PALM due to different grid resolution
    smois_m = np.median(smois_lvl[smois_lvl<1])
    smois_lvl[smois_lvl==1] = smois_m
    return smois_lvl

for d in range(zs.shape[0]):
    smois_land[d,:,:] = calc_soil_moisture(smois[0,d,:,:])
    
for iy in range(smois.shape[2]):
    for ix in range(smois.shape[3]):
        init_soil_t[:,iy,ix] = np.interp(dz_soil, zs, tslb[0,:,iy,ix])
        init_soil_m[:,iy,ix] = np.interp(dz_soil, zs, smois_land[:,iy,ix])


init_soil_tyx = np.empty((dz_soil.shape[0],y.shape[0],x.shape[0]))
init_soil_myx = np.empty((dz_soil.shape[0],y.shape[0],x.shape[0]))

for i in range(dz_soil.shape[0]):
    init_soil_tyx[i,:,:] = rgc.grid_points_change(init_soil_t[i,:,:], x.shape[0], y.shape[0], interp_mode)
    init_soil_myx[i,:,:] = rgc.grid_points_change(init_soil_m[i,:,:], x.shape[0], y.shape[0], interp_mode)




##############################################################################
# Write to NetCDF file
# Based on INIFOR format
print('Writing NetCDF file',flush=True)
nc_output = xr.Dataset()
res_origin = str(dx) + 'x' + str(dy) + ' m'
nc_output.attrs['description'] = 'Contains dynamic data from WRF mesoscale'
nc_output.attrs['history'] = 'Created at ' + time.ctime(time.time())
nc_output.attrs['source']= 'netCDF4 python'
nc_output.attrs['origin_lat'] = np.float(lat_palm)
nc_output.attrs['origin_lon'] = np.float(lon_palm)
nc_output.attrs['z'] = np.float(0)
nc_output.attrs['x'] = np.float(0)
nc_output.attrs['y'] = np.float(0)
nc_output.attrs['rotation_angle'] = np.float(0)
nc_output.attrs['origin_time'] =  str(times[0]) + ' UTC'
nc_output.attrs['end_time'] =  str(times[-1]) + ' UTC'


nc_output['x'] = xr.DataArray(x, dims=['x'], attrs={'units':'m'})
nc_output['y'] = xr.DataArray(y, dims=['y'], attrs={'units':'m'})
nc_output['z'] = xr.DataArray(z, dims=['z'], attrs={'units':'m'})
nc_output['zsoil'] = xr.DataArray(dz_soil, dims=['zsoil'], attrs={'units':'m'})
nc_output['xu'] = xr.DataArray(xu, dims=['xu'], attrs={'units':'m'})
nc_output['yv'] = xr.DataArray(yv, dims=['yv'], attrs={'units':'m'})
nc_output['zw'] = xr.DataArray(zw, dims=['zw'], attrs={'units':'m'})
nc_output['time'] = xr.DataArray(times_sec, dims=['time'], attrs={'units':'seconds'})


nc_output.to_netcdf(f'dynamic_files/{case_name}_dynamic_{ts}')



nc_output['init_soil_m'] = xr.DataArray(init_soil_myx, dims=['zsoil','y','x'], 
         attrs={'units':'m^3/m^3','lod':np.int32(2), 'source':'WRF', 'long_name':'volumetric soil moisture (m^3/m^3)'}) 
nc_output['init_soil_t'] = xr.DataArray(init_soil_tyx, dims=['zsoil','y','x'], 
         attrs={'units':'K', 'lod':np.int32(2), 'source':'WRF', 'long_name':'soil temperature (K)'}) 

# output boundary conditions to PALM input
# directions: 0 left, 1 right
#             0 south, 1 north

nc_output['init_atmosphere_pt'] = xr.DataArray(pt_init,dims=['z'],
         attrs={'units':'K', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_pt'] = xr.DataArray(pt_tmp[:,:,:,0],dims=['time', 'z', 'y'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_pt'] = xr.DataArray(pt_tmp[:,:,:,-1],dims=['time', 'z', 'y'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_pt'] = xr.DataArray(pt_tmp[:,:,0,:],dims=['time', 'z', 'x'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_pt'] = xr.DataArray(pt_tmp[:,:,-1,:],dims=['time', 'z', 'x'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_pt'] = xr.DataArray(pt_tmp[:,-1,:,:],dims=['time', 'y', 'x'],
         attrs={'units':'K', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_qv'] = xr.DataArray(qv_init,dims=['z'],
         attrs={'units':'kg/kg', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_qv'] = xr.DataArray(qv_tmp[:,:,:,0],dims=['time', 'z', 'y'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_qv'] = xr.DataArray(qv_tmp[:,:,:,-1],dims=['time', 'z', 'y'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_qv'] = xr.DataArray(qv_tmp[:,:,0,:],dims=['time', 'z', 'x'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_qv'] = xr.DataArray(qv_tmp[:,:,-1,:],dims=['time', 'z', 'x'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_qv'] = xr.DataArray(qv_tmp[:,-1,:,:],dims=['time', 'y', 'x'],
         attrs={'units':'kg/kg', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_u'] = xr.DataArray(u_init,dims=['z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_u'] = xr.DataArray(u_tmp[:,:,:,0],dims=['time', 'z', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_u'] = xr.DataArray(u_tmp[:,:,:,-1],dims=['time', 'z', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_u'] = xr.DataArray(u_tmp[:,:,0,:],dims=['time', 'z', 'xu'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_u'] = xr.DataArray(u_tmp[:,:,-1,:],dims=['time', 'z', 'xu'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_u'] = xr.DataArray(u_tmp[:,-1,:,:],dims=['time', 'y', 'xu'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_v'] = xr.DataArray(v_init,dims=['z'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_v'] = xr.DataArray(v_tmp[:,:,:,0],dims=['time', 'z', 'yv'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_v'] = xr.DataArray(v_tmp[:,:,:,-1],dims=['time', 'z', 'yv'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_v'] = xr.DataArray(v_tmp[:,:,0,:],dims=['time', 'z', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_v'] = xr.DataArray(v_tmp[:,:,-1,:],dims=['time', 'z', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_v'] = xr.DataArray(v_tmp[:,-1,:,:],dims=['time', 'yv', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})

nc_output['init_atmosphere_w'] = xr.DataArray(w_init,dims=['zw'],
         attrs={'units':'m/s', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_left_w'] = xr.DataArray(w_tmp[:,:,:,0],dims=['time', 'zw', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_right_w'] = xr.DataArray(w_tmp[:,:,:,-1],dims=['time', 'zw', 'y'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_south_w'] = xr.DataArray(w_tmp[:,:,0,:],dims=['time', 'zw', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_north_w'] = xr.DataArray(w_tmp[:,:,-1,:],dims=['time', 'zw', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_top_w'] = xr.DataArray(w_tmp[:,-1,:,:],dims=['time', 'y', 'x'],
         attrs={'units':'m/s', 'source':'WRF', 'res_origin':res_origin})

nc_output['surface_forcing_surface_pressure'] = xr.DataArray(pres_init,dims=['time'],
         attrs={'units':'Pa', 'lod':np.int32(1), 'source':'WRF', 'res_origin':res_origin})


nc_output['ls_forcing_ug'] = xr.DataArray(geo_wind_u_f,dims=['time','z'],
         attrs={'units':'m/s', 'long_name':'u wind component geostrophic', 'source':'WRF', 'res_origin':res_origin})
nc_output['ls_forcing_vg'] = xr.DataArray(geo_wind_v_f,dims=['time','z'],
         attrs={'units':'m/s', 'long_name':'v wind component geostrophic', 'source':'WRF', 'res_origin':res_origin})



for var in nc_output.data_vars:
    encoding = {var: {'dtype': 'float32', '_FillValue': -9999, 'zlib':True}}
    nc_output[var].to_netcdf(f'dynamic_files/{case_name}_dynamic_{ts}', encoding=encoding, mode='a')


print('Add to your *_p3d file the: ' + '\n soil_temperature = ' + repr(init_soil_t.mean(axis=1).mean(axis=1)) +
      '\n soil_moisture = ' + repr(init_soil_m.mean(axis=1).mean(axis=1)) + '\n deep_soil_temperature = ' + repr(init_soil_tmn))

end = datetime.now()
print('PALM dynamic input file is ready. Duration: {}'.format(end - start))
print('Start time: '+str(times[0]))
print('End time: '+str(times[-1]))
print('Time step: '+str(time_step_sec)+' seconds')
del u_int, v_int, w_int, qv_int, pt_int, pres_int, tk_int
