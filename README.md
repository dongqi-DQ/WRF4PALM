# WRF4PALM v1.1 

*if you wish to use WRF4PALM v1.0, please go to https://github.com/dongqi-DQ/WRF4PALM/tree/v1.0*

## Contents
1. [What's new?](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1/#whats-new-in-v11)
2. [Instructions](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1/#instrustions)
3. [Quick comparison script](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1/#quick-compare-wrf--palm)

## what's new in v1.1?
- move from wrf-python to salem to modify RAM usage and computation time 
- use xarray instead of netCDF4 package to modify RAM usage and computation time
- apply multiprocessing to improve computation time 
- users now only need to edit namelist instead of editing the script  
- add surface variables (e.g. U10, V10, T2, and Q2) for surface NaN solver 
- read WRF projection info when locate PALM domain 
- allow users to specify the projection of PALM simulation
- geostrophic winds are estimated using geopotential height instead of pressure 

## Instrustions
**How to use WRF4PALM v1.1##
1. Download the entire source code to local
2. Provide your own WRF output
3. Edit [namelist](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#namelist)
4. [Run WRF4PALM](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#one_line_command)

### namelist 

### One line command

## Quick compare WRF & PALM


# Remark  
water temperature in static driver
we might release a static driver generator using global data set from Google earth engine

We encourage WRF4PALM users to use the GitHub **Issue** system if they encountered issues or problems using WRF4PALM.
