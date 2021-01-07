# WRF4PALM v1.0.1   

**Tools to create a dyanmic driver to enable WRF-PALM offline nesting**
------
**[Update info v1.0.1]**

To make the process of WRF4PALM more straight forward for users, some updates have been made in v1.0.1:  
1. Users now only need to edit `namelist.dynamic` instead of running two different Python scripts
2. The two Python scripts `create_cfg.py` and `create_dynamic.py` are combined into one script to save manual work
3. Modify geostrophic winds calculation. 

------

This repository contains the following files:

- `namelist.dynamic` - Python namelist for user input 
- `create_dynamic.py` - the major WRF4PALM script which interpolates WRF output to a PALM dynamic driver   
- scripts in `util` - functions to use in WRF4PALM  
- Example cfg files in `cfg_input`  
- Example dynamic drivers in `dynamic_files` 

## Note:
- The scripts may only work on Python 3. Some features or functions may need some tweaks if the scripts are used in Python2. 
- When processing to high resolution and large PALM domain, the script may cost a lot of RAM.   
- The `interplevel` function from the [wrf-python](https://wrf-python.readthedocs.io/en/latest/) package may have bugs (creating NaN) when interpolate to certain height levels. Several functions have been applied in the coupler to avoid the NaN values, which may slow down the surface NaN solver.
- To minimise interpolation errors, the WRF projection should be centred close to the centre of the PALM domain.

# Step 1: Edit namelist.dynamic

## Tips if you are new to Python namelist:  
1. Square brakets indicate the name/key of each section, such as `[case]`, `[domain]`, and `[stretch]`. Please do not change these keys.  
2. Under each section, there are parameters for users to specify. Please include commas (,) at the end of each line. Otherwise, WRF4PALM may not work.  
3. Use double quotes for strings.  
4. Comments in the namelist still begin with a hash mark (#). 

## case
Users first need to specify the case name for the dynamic driver in `[case]` section:

```python
[case]
case_name = "your case name"     # case name as you prefer
```

## domain
Specify PALM domain configuration under `[domain]` section:

```python
[domain]
centlat   = -43.529469,                           # latitude of domain centre
centlon   = 172.599316,                           # longitude of domain centre
nx        = 200,                                  # number of grid points along x-axis
ny        = 200,                                  # number of grid points along y-axis
nz        = 200,                                  # number of grid points along z-axis
dx        = 80.0,                                 # number of grid points along x-axis
dy        = 80.0,                                 # number of grid points along y-axis
dz        = 16.0,                                 # number of grid points along z-axis
z_origin  = 0.0,                                  # elevated mean grid position (elevated terrain)
```

Here `z_origin` can be specified to consider elevated mean grid position (elevated terrain) thus avoiding many NAN and unrealistic extrapolation towards z=0 m.  
If `z_origin` is not known, leave it as 0 m. PALM will automatically exclude all data below the elevated level if terrain is included in the simulation. 

## strech
The `[strech]` section is to enable a vertically streched grid. The parameters are identical to PALM's parameters. 
```python
[stretch]
dz_stretch_factor = 1.0,        # stretch factor for a vertically stretched grid
                                # set this to 1.0 if no streching required
dz_stretch_level = 1200.0,      # Height level above which the grid is to be stretched vertically (in m)

dz_max = 30.0,                  # allowed maximum vertical grid spacing (in m)
```

## wrf
The `[wrf]` section includes all the parameters related to WRF output. Users must specify:

1. WRF output flie name and location
```python
wrf_output = "wrf_output/your_wrf_output.nc",
```

2. Interpolation method (`"linear"` or `"nearest"`)
```python
interp_mode = "linear",
```

3. Specify the start and end time stamp as well as the update frequency of boundary conditions to be used in PALM simulation
```python
start_year = 2019,
start_month = 7,
start_day = 1,
start_hour = 0,

end_year = 2019,
end_month = 7,
end_day = 1,
end_hour = 10,

interval = 6,
ts = '1hour',
```


**Remark:**   
- here users can define the start/end time and update frequency depending on their own requirments. For example, a user has hourly WRF output, if a user only wants the boundary conditions update every 2 hours. The `interval` is set to 2 and `ts` is 2hour, which gives a dynamic driver with filename `case_name_2hour`. Here I have 10-minute WRF output and want hourly update in PALM, then `interval = 6` and `ts = '1hour'` (simple math :).  

- When the update frequency is not a divisor of the total run time, the script will round up the last step of time to avoid automatic termination of PALM. For example, if one has a 24-hour PALM simulation (`end_time = 86400`) with 5-hour update of boundary conditions, the time stamp in the dynamic file should be 0s, 18000s, 36000s, 54000s, 72000s, 86400s (instead of 90000s because 86400s may reach the final time step in WRF output). If the 86400s step is not in the dynamic file, PALM may return some errors and terminate at the 72000s step. 

## soil
In `[soil]` section, users must specify the depth of soil layers. The default settings are identical to PALM's default 8 layers.
```python
[soil]
# layers for soil temperature and moisture calculation
# this shall be changed depending on different cases

dz_soil = 0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86,
```

# Step 2: Run create_dynamic.py

Now run the script `create_dynamic.py`. If it is successfully executed, the dynamic file will be stored in `dynamic_files`. 

The Python script will return progress information while running:
```bash
Loading WRF netCDF: wrf_output/wrfout_d04_2017-02-10_00.nc
WRF output reading done.
Interpolating horizontal fields...
100%|##########| 7/7 [00:55<00:00,  7.96s/it]Done.

Interpolating unstaggered vertical levels: 200it [10:04,  3.02s/it]
Interpolating staggered vertical levels: 199it [02:03,  1.62it/s]
Vertical interpolation done.

Resolving surface NaNs: 100%|##########| 7/7 [01:49<00:00, 15.70s/it]

Calculating soil temperature and moisture from WRF
Writing NetCDF file
```

It will give the soli moisture and soil temperature information for PALM simulation after the dynamic driver is created:
```bash
Add to your *_p3d file the: 
 soil_temperature = array([288.28391647, 288.28391647, 288.28391647, 288.38816098,
       289.22211713, 290.37291891, 290.48805698, 289.21989822])
 soil_moisture = array([0.18927199, 0.18927199, 0.18927199, 0.18951527, 0.19146148,
       0.19428751, 0.19848686, 0.23685172])
 deep_soil_temperature = 285.82938
PALM dynamic input file is ready. Script duration: 0:15:07.212707
Start time: 2019-07-01 00:00:00
End time: 2019-07-01 10:00:00
Time step: 3600.0 seconds
```
The information above shows that it took 15 minutes to finish processing. The start time, end time and time step are also given showing data in the dynamic driver is as desired. If all the information is correct, then the dynamic driver is ready to use for PALM.

------ End of README ------

Development of WRF4PALM is based on WRF2PALM (https://github.com/ricardo88faria/WRF2PALM). 

A full documentation is still under construction, if you have any queries please contact the author.

**Contact: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)**
