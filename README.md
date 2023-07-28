# WRF4PALM v1.1  [![DOI](https://zenodo.org/badge/258736274.svg)](https://zenodo.org/badge/latestdoi/258736274)


*If you wish to use WRF4PALM v1.0, please go to [WRF4PALM v1.0](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.0)*

**Want to create your own static driver? Check out [GEO4PALM](https://github.com/dongqi-DQ/GEO4PALM)!**

## Contents
1. [What's new?](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1/#whats-new-in-v11)
2. [Instructions](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1/#instrustions)
3. [Quick comparison script](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1/#quick-compare-wrf--palm)

## What's new in v1.1.2?
1. Users now have option to calculate geostropihc wind using geopotential height (option `"z"`) or pressure (option `"p"`) via `geostrophic` option in `[case]`.  
2. The surface NaN solver now is taking surface variables from the correct altitude by including terrain height in WRF.

## What's new in v1.1?
- multiple WRF output files are allowed
- move from wrf-python to salem to modify RAM usage and computation time
- use xarray instead of netCDF4 package to modify RAM usage and computation time
- apply multiprocessing to improve computation time
- users now only need to edit namelist instead of editing the script  
- add surface variables (e.g. U10, V10, T2, and Q2) for surface NaN solver
- read WRF projection info when locate PALM domain
- allow users to specify the projection of PALM simulation
- geostrophic winds are estimated using geopotential height instead of pressure

## Instrustions
**How to use WRF4PALM v1.1**
1. Download the entire source code to local
2. Provide your own WRF output
3. Edit [namelist](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#namelist)
4. [Run WRF4PALM](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#one-line-command)


### namelist
In v1.1, users don't have to edit the main script, and only need to edit the namelist file to provide their input (for examples please see `namelist.wrf4palm`).

There are 5 sections in the namelist:
- [case](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#case): users to provide case name and multiprocessing information
- [domain](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#domain): PALM domain configuration
- [stretch](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#stretch): if a vertically streched grid is used
- [wrf](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#wrf): information about WRF output and start/end time in the dynamic driver
- [soil](https://github.com/dongqi-DQ/WRF4PALM/tree/v1.1#soil): soil layers and dummy soil moisture information

#### case
In the `case` section, users need to provide their case name and the maximum number of CPUs they want to use in WRF4PALM (here the number is 4).
```
[case]
case_name = "wrf4palm_test", # specify your case name here
max_pool = 4,                # specify the maximum number of CPUs to use
```

#### domain
In the `domain` section, users need to provide PALM domain configuration (dx, dy, dz, nx, ny, nz, and z_origin), the latitude and longitude at PALM domain centre, and the projection of PALM domain. The projection of PALM domain and centre lat/lon are used to locate PALM domain in the WRF domain. The projection of PALM domain should be identical to the projection of PALM static driver, if the user has one. If users do not have the projection information, they can leave the field empty as `palm_proj = "",` such that WRF4PALM v1.1 will use the projetion of WRF directly.

```
[domain]
palm_proj = "EPSG:2193",    # projection of PALM
centlat   = -35.7853,       # latitude of domain centre
centlon   = 174.1,          # longitude of domain centre
nx        = 200,            # number of grid points along x-axis
ny        = 200,            # number of grid points along y-axis
nz        = 120,            # number of grid points along z-axis
dx        = 50.0,           # number of grid points along x-axis
dy        = 50.0,           # number of grid points along y-axis
dz        = 10.0,           # number of grid points along z-axis
z_origin  = 0.0,            # elevated mean grid position (elevated terrain)
```

#### stretch
In the `stretch` section, useres can define vertically streched grid spacing. The parameters are identical to those in PALM. If no streching is required, leave `dz_stretch_factor=1.0,`
```
[stretch]
dz_stretch_factor = 1.0,        # stretch factor for a vertically stretched grid
                                # set this to 1.0 if no streching required
dz_stretch_level = 1200.0,      # Height level above which the grid is to be stretched vertically (in m)

dz_max = 30.0,                  # allowed maximum vertical grid spacing (in m)
```

#### wrf
WRF4PALM users must provide their own WRF output. Users must specify the directory (`wrf_path`) to access WRF netcdf output files, and WRF output filenames. WRF4PALM v1.1 allows users to provide one or multiple WRF files. Users can either provide a list of filenames, e.g.:  
`wrf_output = "wrfout_d04_2020-12-25_12-00-00", "wrfout_d04_2020-12-26_12-00-00"`
or a string glob in the form:
`wrf_output = "wrfout_d04_2020-12-*", `

Users also need to specify the interpolation mode (`interp_mode`) to interpolate WRF output onto PALM grid. Both `"linear"` and `"nearest"` are allowed, while we recommend using `"linear"`.

The start and end datetime of PALM simulation must be provided. The PALM dynamic driver update frequency is controlled by `dynamic_ts` (unit: seconds), e.g. `dynamic_ts = 3600.0,` means the boudnary conditions will be updated every hour.

```
[wrf]
wrf_path = "./wrf_output/",
wrf_output = "wrfout_d04_2020-12-25_12-00-00", "wrfout_d04_2020-12-26_12-00-00",

interp_mode = "linear",

start_year = 2020,
start_month = 12,
start_day = 25,
start_hour = 13,

end_year = 2020,
end_month = 12,
end_day = 26,
end_hour = 10,

dynamic_ts = 3600.0,         # PALM dynamic driver update frequency (seconds)

```

**Note**: leading zeros are not permitted in the datetime configuration. For example, if the `start_month` is January, then the namelist should have `start_month = 1,` instead of `start_month = 01,`.

#### soil
In the `soil` section, users need to config the soil layers (`dz_soil`). In case when soil moisture output in WRF is all zeros (due to WRF's parameterisation), a dummy value can be chosen (e.g. `msoil = 0.3,`).

```
[soil]
# layers for soil temperature and moisture calculation
# this shall be changed depending on different cases

dz_soil = 0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86,
msoil = 0.3,         # dummy value in case soil moisture from WRF output is 0.0
```

### One line command
Once the namelist is ready, users can run WRF4PALM using the one line command:
```
python run_config_wrf4palm.py [your namelist]
```

**Execution example**
```
Reading WRF
cfg file is saved: wrf4palm_test
Start horizontal interpolation
Calculating soil temperature and moisture from WRF
100%|█████████████████████████████████████████████████████████████████| 200/200 [00:29<00:00,  6.74it/s]
Start vertical interpolation
Processing QVAPOR for west and east boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:23<00:00,  5.13it/s]
Processing QVAPOR for south and north boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:23<00:00,  5.14it/s]
Processing pt for west and east boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:25<00:00,  4.68it/s]
Processing pt for south and north boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:25<00:00,  4.73it/s]
Processing W for west and east boundaries
100%|█████████████████████████████████████████████████████████████████| 119/119 [00:23<00:00,  5.00it/s]
Processing W for south and north boundaries
100%|█████████████████████████████████████████████████████████████████| 119/119 [00:23<00:00,  5.08it/s]
Processing U for west and east boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:22<00:00,  5.26it/s]
Processing U for south and north boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:23<00:00,  5.19it/s]
Processing V for west and east boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:23<00:00,  5.19it/s]
Processing V for south and north boundaries
100%|█████████████████████████████████████████████████████████████████| 120/120 [00:23<00:00,  5.20it/s]
Processing top boundary conditions...
100%|█████████████████████████████████████████████████████████████████████| 5/5 [01:08<00:00, 13.70s/it]
Geostrophic wind estimation...
100%|███████████████████████████████████████████████████████████████████| 22/22 [00:02<00:00,  8.80it/s]
Resolving surface NaNs...
100%|█████████████████████████████████████████████████████████████████████| 5/5 [01:05<00:00, 13.14s/it]
Writing NetCDF file
Add to your *_p3d file:
 soil_temperature = [288.09888275146477, 288.42630248766665, 289.2634380215685, 289.9926226307226, 292.3663809474804, 293.1886674499512, 293.1886674499512, 293.1886674499512]
 soil_moisture = [0.29, 0.29, 0.29, 0.29, 0.29, 0.29, 0.29, 0.29]
 deep_soil_temperature = 294.6415

PALM dynamic input file is ready. Script duration: 0:07:26.348010
Start time: 2020-12-25T13:00:00.000000000
End time: 2020-12-26T10:00:00.000000000
Time step: 3600.0 seconds
```

If the execution is successful, the dynamic driver will be ready in `dynamic_files` with the `case_name` and start timestamp user specified. A cfg reference file will also be stored in `cfg_files` which contains domain configuration and soil temperatuer and moisture information. An example dynamic driver and an example cfg file are provided in `dynamic_files` and `cfg_files`, respectively.

## Quick compare WRF & PALM

In order for users to quickly check the quality of the dynamic driver generated by WRF4PALM, we provide a quick comparison script. Five variables are allowed (can be in uppercase or lowercase):
- U
- V
- W
- PT
- QV

Three plot types are provided:
1. **zcross**: vertical cross sections of west/east/south/north boundaries for the user specified variable and timestamp
```
python3 quick_compare.py [your namelist] zcross [variable name]
```
then the script will ask for the timestamp:
```
Please enter the timestamp (yyyy-mm-dd-hh):
```
Once the timestamp is given, the script will return a comparison plot.

2. **pr**: vertical profiless of west/east/south/north boundaries for the user specified variable and timestamp
```
python3 quick_compare.py [your namelist] pr [variable name]
```
then the script will ask for the timestamp:
```
Please enter the timestamp (yyyy-mm-dd-hh):
```
Once the timestamp is given, the script will return a comparison plot.  
Note that the vertical profiles are horizontally averaged and hence the comparison only gives a approximate reference regarding the performance of WRF4PALM.

2. **ts**: time series of west/east/south/north boundaries for the user specified variable and altitude
```
python3 quick_compare.py [your namelist] ts [variable name]
```
then the script will ask for the altitude:
```
Please enter the vertical level in m:
```
Once the vertical level is given, the script will return a comparison plot.  
Note that the time series are horizontally averaged and hence the comparison only gives a approximate reference regarding the performance of WRF4PALM.

## Remark
- [`Surface_NaN_Solver.pdf`](https://github.com/dongqi-DQ/WRF4PALM/blob/v1.1/Surface_NaN_Solver.pdf) provides a short documentation explaining how the surface nans are resolved.
- The WRF4PALM v1.1 python environemnt is available in [`wrf4palm_env.yml`](https://github.com/dongqi-DQ/WRF4PALM/blob/v1.1/wrf4palm_env.yml).

# Note  
- We noticed that PALM uses a water temperature of 283 K as default, which may lead to a stable layer over water bodies (if there are any in the PALM simulation domain). We recommend users to modify the water temperatuer using the static driver.
- We may release a static driver generator using global data set from Google earth engine and SST from ERA5 (date TBC).
- Geostrophic winds are only an estimation while the accuracy of the estimation still needs further discussion and investigation. This problem is the same in INIFOR.
- We encourage WRF4PALM users to use the GitHub **Issue** system if they encountered any issues or problems using WRF4PALM such that communications and trouble shooting will be easier.

--------------------------------------------------------------------------------------------
### End of README
--------------------------------------------------------------------------------------------

Development of WRF4PALM is based on WRF2PALM (https://github.com/ricardo88faria/WRF2PALM).

A full documentation is still under construction, if you have any queries please contact the author or open a new issue.

--------------------------------------------------------------------------------------------
**Contact: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)**

**How to cite**

Lin, D., Khan, B., Katurji, M., Bird, L., Faria, R., and Revell, L. E.: WRF4PALM v1.0: a mesoscale dynamical driver for the microscale PALM model system 6.0, Geosci. Model Dev., 14, 2503–2524, https://doi.org/10.5194/gmd-14-2503-2021, 2021.
