#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------#
# WRF4PALM
#
# functions to resolve NaN values near surface
#
# @author: Dongqi Lin (dongqi.lin@canterbury.ac.nz)
#--------------------------------------------------------------------------------#
import numpy as np

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
                    ds_we[var].loc[dict(time=ts)][:,j,bc] = surface_nan_w(ds_we[var].sel(time=ts)[:,j,bc].data, ds_we.zw.data)
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
                    ds_sn[var].loc[dict(time=ts)][:,bc,i] = surface_nan_w(ds_sn[var].sel(time=ts)[:,bc,i].data, ds_sn.zw.data)
                else:
                    surface_var = surface_var_dict[var]
                    surface_value = surface_var.sel(time=ts)[bc, i].data
                    ds_sn[var].loc[dict(time=ts)][:,bc,i] = surface_nan_s(ds_sn[var].sel(time=ts)[:,bc,i].data, z, surface_value)
    return var, (ds_we[var], ds_sn[var])

Z0_SURFACE = 0.5   # roughness length (m) for the near-surface wind profile: WRF's urban value


def _ground(zlev):
    """Height of the PALM ground from scalar levels z (z[0] is half a spacing above it)."""
    return zlev[0] - 0.5 * (zlev[1] - zlev[0])


def surface_nan_s(data, z, s2):
    """Fill the NaN levels at the bottom of a scalar column (PALM levels below the lowest
    valid WRF level): a straight line in height through the WRF 2 m value s2, placed 2 m
    above the PALM ground, and the first valid level. The level below 2 m gets s2 itself.
    (Patched: the original used a slope of (s2 - s_first) / 2 anchored at the first valid
    level, which is far too steep when that level is 20 m up.)"""
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0:
        k1 = int(nan_idx.max()) + 1
        if k1 >= data.size:
            return data
        z2 = _ground(z) + 2.0
        a = (data[k1] - s2) / (z[k1] - z2) if z[k1] > z2 else 0.0
        for i in range(k1):
            data[i] = s2 if z[i] <= z2 else s2 + a * (z[i] - z2)
    return data


def surface_nan_uv(data, z, uv10):
    """Fill the NaN levels at the bottom of a wind-component column: below 10 m a logarithmic
    profile anchored at the WRF 10 m wind uv10 (roughness length Z0_SURFACE), between 10 m
    and the first valid level a straight line. If the first valid level is itself below 10 m,
    the logarithmic profile is anchored at that level. The wind vanishes at the ground.
    (Patched: the original fitted a logarithm through the 10 m value and the first valid
    level as if they were 10 m apart.)"""
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0:
        k1 = int(nan_idx.max()) + 1
        if k1 >= data.size:
            return data
        ground = _ground(z)
        h1 = z[k1] - ground
        if h1 > 10.0:
            ha, ua = 10.0, uv10
        else:
            ha, ua = h1, data[k1]
        for i in range(k1):
            h = z[i] - ground
            if h <= ha:
                data[i] = 0.0 if h <= Z0_SURFACE else ua * np.log(h / Z0_SURFACE) / np.log(ha / Z0_SURFACE)
            else:
                data[i] = ua + (data[k1] - ua) * (h - ha) / (h1 - ha)
    return data


def surface_nan_w(data, zw=None):
    """Fill the NaN levels at the bottom of a vertical-velocity column: w goes linearly to
    zero at the ground from the first valid level. Without zw (the w levels), a constant
    fill with the first valid value. (Patched: the original copied the value two levels up.)"""
    nan_idx = np.argwhere(np.isnan(data))
    if nan_idx.size > 0:
        k1 = int(nan_idx.max()) + 1
        if k1 >= data.size:
            return data
        if zw is None or len(zw) < 2:
            data[:k1] = data[k1]
        else:
            ground = zw[0] - (zw[1] - zw[0])
            for i in range(k1):
                data[i] = data[k1] * (zw[i] - ground) / (zw[k1] - ground)
    return data
