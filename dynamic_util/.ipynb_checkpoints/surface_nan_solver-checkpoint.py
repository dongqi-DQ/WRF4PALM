#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# WRF4PALM
#
# functions to resolve NaN values near surface
#
# @author: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)
#--------------------------------------------------------------------------------#
import numpy as np

def surface_nan_uv(data,z, uv10):
    '''
    assuming U = a*log(z)+b
    '''
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0 :
        nan_idx = int(nan_idx[-1])+1
        a = (data[nan_idx]-uv10)/np.log(z[nan_idx]-10)
        b = uv10-a*np.log(10)
        for i in range(0,int(nan_idx)):
            data[i] = a*np.log(z[i])+b
    return(data)
def surface_nan_s(data,z,s2):
    '''
    reslove surface nan for scalars
    assumming s = a*z+b
    '''
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0 :
        nan_idx = int(nan_idx[-1] + 1)
        a = (data[nan_idx]- s2)/(z[nan_idx]-2)
        b = s2-a*2
        for i in range(0,int(nan_idx)):
            data[i] = a*z[i]+b
    return(data)
def surface_nan_w(data):
    '''
    do nothing for vertical wind 
    '''
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0 :
        nan_idx = int(nan_idx[-1] + 1)
        data[:nan_idx] = data[nan_idx+1]
    return(data)

def solve_surface(all_ts, ds_we, ds_sn, surface_var_dict, var):
    z = ds_we.z.data
    for ts in all_ts:
        for j in range(0,ds_we[var].shape[2]):
            for bc in [0,-1]:
                if var=="U" or var=="V":
                    surface_var = surface_var_dict[var]
                    surface_value = surface_var.sel(time=ts)[j,bc].data
                    ds_we[var].loc[dict(time=ts)][:,j,bc] = surface_nan_uv(ds_we[var].sel(time=ts)[:,j,bc].data,z,surface_value)
                elif var=="W":
                    ds_we[var].loc[dict(time=ts)][:,j,bc] = surface_nan_w(ds_we[var].sel(time=ts)[:,j,bc].data)
                else:
                    surface_var = surface_var_dict[var]
                    surface_value = surface_var.sel(time=ts)[j,bc].data
                    ds_we[var].loc[dict(time=ts)][:,j,bc] = surface_nan_s(ds_we[var].sel(time=ts)[:,j,bc].data, z, surface_value)

        for i in range(0,ds_sn[var].shape[3]):
            for bc in [0,-1]:
                if var=="U" or var=="V":
                    surface_var = surface_var_dict[var]
                    surface_value = surface_var.sel(time=ts)[bc, i].data
                    ds_sn[var].loc[dict(time=ts)][:,bc,i] = surface_nan_uv(ds_sn[var].sel(time=ts)[:,bc,i].data,z,surface_value)
                elif var=="W":
                    ds_sn[var].loc[dict(time=ts)][:,bc,i] = surface_nan_w(ds_sn[var].sel(time=ts)[:,bc,i].data)
                else:
                    surface_var = surface_var_dict[var]
                    surface_value = surface_var.sel(time=ts)[bc, i].data
                    ds_sn[var].loc[dict(time=ts)][:,bc,i] = surface_nan_s(ds_sn[var].sel(time=ts)[:,bc,i].data, z, surface_value)
    return var, (ds_we[var], ds_sn[var])

