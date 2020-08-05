# WRF4PALM 
**Tools to create a dyanmic driver to enable WRF-PALM offline nesting**


This repository contains the following files:

- `create_cfg.py` - Python script to create a cfg file including PALM domain configuration 
- `create_dynamic.py` - the major WRF4PALM script which interpolates WRF output to a PALM dynamic driver   
- scripts in `util` - functions to use in WRF4PALM  
- Example cfg files in `cfg_input`  
- Example dynamic drivers in `dynamic_files` 

## Note:
- The scripts may only work on Python 3. Some features or functions may need some tweaks if the scripts are used in Python2. 
- When coupling to high resolution and large PALM domain, the coupling script may cost a lot of RAM.   
- The `interplevel` function from the [wrf-python](https://wrf-python.readthedocs.io/en/latest/) package may have bugs (creating NaN) when interpolate to certain height levels. Several functions have been applied in the coupler to avoid the NaN values, which may slow down the surface NaN solver.

# Let's get start
## Step 1 Locate the Domain

The users first need to give the domain information in the `create_cfg.py` script. The information including:  

```python
case_name_d01 = 'chch_NW_10m'     # case name as you prefer, but should be consistent with the one used in dynamic script
centlat_d01   = -43.487           # latitude of domain centre
centlon_d01   = 172.537           # longitude of domain centre
dx_d01        = 10                # grid spacing in meters along x-axis
dy_d01        = 10                # grid spacing in meters along y-axis
dz_d01        = 10                # grid spacing in meters along z-axis
nx_d01        = 360               # number of grid points along x-axis
ny_d01        = 360               # number of grid points along y-axis
nz_d01        = 360               # number of grid points along z-axis
```

Run the script `create_cfg.py` and then you will get a cfg file (in `cfg_input` folder) containing the raw domain configuration for `create_dynamic.py`.

## Step 2 Process WRF for PALM

1. Specify case name, which should be the same as the one specified in Step 1. 
```python
case_name = 'chch_NW_10m' # case name as you specified in create_cfg.py
```

2. Specify the WRF output file (in `wrf_output` folder) to process 
```python
wrf_file = 'wrf_output/your_wrf_output_file_name' 
```

3. Specify the start and end time stamp as well as the update frequency of boundary conditions you want to used in PALM simulation  

```python
dt_start    = datetime(2017, 2, 11, 20,)  #  start time in YYYY/MM/DD/HH format
dt_end      = datetime(2017, 2, 12, 2,)   #  end time in YYYY/MM/DD/HH format
interval    = 1                           #  define time interval of WRF output to be read for the coupling
ts          = '1hour'                     #  specify the update frequency of boundary conditions which will 
                                          #  show in the dynamic input filename
                                          #  this works as a reference in case the update frequency calculation went wrong
```

**Remark:**   
- here users can define the start/end time and update frequency depending on their own requirments. For example, I have hourly WRF output, if I only want the boundary conditions update every 2 hours. The `interval` is set to 2 and `ts` is 2hour, which gives a dynamic driver with filename `case_name_2hour`. If you have 10-minute WRF output and want hourly update in PALM, then `interval = 6` and `ts = '2hour'` (simple math :).  
- When the update frequency is not a divisor of the total run time, the script will round up the last step of time to avoid automatic termination of PALM. For example, if one has a 24-hour PALM simulation (`end_time = 86400`) with 5-hour update of boundary conditions, the time stamp in the dynamic file should be 0s, 18000s, 36000s, 54000s, 72000s, 86400s (instead of 90000s because 86400s may reach the final time step in WRF output). If the 86400s step is not in the dynamic file, PALM may return some errors and terminate at the 72000s step. 


4. Specify the depth of 8 soil layers 

```python
dz_soil = np.array([0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86]) # this is the default setup in PALM
```

5. Specify `z_origin` to consider elevated mean grid position (elevated terrain) thus avoiding many NAN and unrealistic extrapolation towards z=0. 
```python
z_origin = 0
```
If `z_origin` is not known, leave it as 0 m. PALM will automatically exclude all data below the elevated level if terrain is included in the simulation.

6. (**Optional**) Define vertically streched grid spacing. The parameters are identical to those in PALM.
```python
dz_stretch_factor = 1.02    # Stretch factor for a vertically stretched grid. Set to 1 if no strech required.

dz_stretch_level = 1200     # Height level above which the grid is to be stretched vertically (in m)

dz_max = 30                 # aloowed maximum vertical grid spacing (in m)
```

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
Start time: 2017-02-11 20:00:00
End time: 2017-02-12 02:00:00
Time step: 3600.0 seconds
```
The information above shows that it took 15 minutes to finish processing. The start time, end time and time step are also given showing data in the dynamic driver is as desired. If all the information is correct, then the dynamic driver is ready to use for PALM.

* * * End of README * * *

Development of WRF4PALM is based on WRF2PALM (https://github.com/ricardo88faria/WRF2PALM).

**Contact: Dongqi Lin (dongqi.lin@pg.canterbury.ac.nz)**
