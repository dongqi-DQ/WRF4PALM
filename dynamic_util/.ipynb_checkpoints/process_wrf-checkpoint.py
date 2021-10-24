#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from tqdm import tqdm
from functools import partial
from multiprocess import Pool
def zinterp( ds, var,lvl ):
    data = ds.salem.wrf_zlevel(var, levels=lvl, use_multiprocessing=False)
    return (lvl, data)

def multi_zinterp(max_pool, ds_in, var, zcoord, ds_out):
    with Pool(max_pool) as p:
        pool_outputs = list(
            tqdm(
                p.imap(partial(zinterp, ds_in, var),zcoord), total=len(zcoord),
                position=0, leave=True
            )
        )
    p.join()
    ## convert dictionary back to dataset
    pool_dict = dict(pool_outputs)
    for lvl in zcoord:
        if var == "W":
            ds_out[var].loc[dict(zw=lvl)] = pool_dict[lvl] 
        else:
            ds_out[var].loc[dict(z=lvl)] = pool_dict[lvl] 
    return ds_out[var]


def process_top(all_ts,top_dict, var):
    ds = top_dict[var][0]
    top_arr = top_dict[var][1]
    zcoord = top_dict[var][2]
    for ts in range(0,len(all_ts)):
        top_arr[ts,:,:] = ds.isel(time=ts).salem.wrf_zlevel(var, levels=zcoord[-1])
    return var, top_arr
