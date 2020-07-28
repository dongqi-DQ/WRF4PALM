# WRF4PALM 
**(for users who don't need to create new static files)**


This repository contains the following files:

- `create_cfg.py` - Python script to create a cfg file including PALM domain information  
- `create_dynamic.py` - the WRF-PALM coupling script  
- scripts in `util` - functions to use in the coupling  
- An example cfg file in `cfg_input`  
- An example dynamic input in `dynamic files` 
- [`WRF-PALM-framewrok.pdf`](https://github.com/dqbuhtig/WRF-PALM-no-static/blob/master/WRF-PALM-framework.pdf) explaining the framework

## Note:
- The scripts may only work on Python 3. Some features or functions may need some tweaks by the Python 2 users themselves. 
- Now users can downscale WRF to PALM when the PALM domain is smaller than one grid spacing of WRF, but this may not be very accurate.  
- Small shifts may apply when locating small domains using latitudes and logitudes.  
- When coupling to high resolution and large PALM domain, the coupling script may cost a lot of RAM. We hence recommend to apply WRF-PALM offline nesting to a large and coarse-resolution PALM domain first and then apply PALM self-nesting for all the small high-resolution domains.  
- The `interplevel` function from the [wrf-python](https://wrf-python.readthedocs.io/en/latest/) package may have bugs (creating NaN) when interpolate to certain height levels. Several functions have been applied in the coupler to avoid the NaN values, but we cannot garantee the same issue wouldn't appear again.

# Let's get start
## Step 1 Locate the Domain

The users first need to give the domain information in the `create_cfg.py` script. The information including [**start from Line 54**](https://github.com/dqbuhtig/WRF-PALM-no-static/blob/master/create_cfg.py#L54):  

```python
case_name_d01 = 'chch_50m_example' # case name as you prefer, but should be consistent with the one used in dynamic script
centlat_d01   = -43.508760         # latitude of domain centre
centlon_d01   = 172.664099         # longitude of domain centre
dx_d01        = 50                 # resolution in meters along x-axis
dy_d01        = 50                 # resolution in meters along y-axis
dz_d01        = 50                 # resolution in meters along z-axis
nx_d01        = 96                 # number of grid points along x-axis
ny_d01        = 96                 # number of grid points along y-axis
nz_d01        = 96                 # number of grid points along z-axis
```

Run the script `create_cfg.py` and then you will get a cfg file (in `cfg_input` folder) containing the domain information for the coupler script `create_dynamic.py`.

## Step 2 Process WRF for PALM

1. Specify case name which should be the same as the one specified in Step 1. [Line 42 in `create_dynamic.py`](https://github.com/dqbuhtig/WRF-PALM-no-static/blob/master/create_dynamic.py#L42)  
```python
case_name = 'chch_50m_example' # case name as you specified in create_cfg.py
```

2. Specify the WRF output file to process [Line 44 in `create_dynamic.py`](https://github.com/dqbuhtig/WRF-PALM-no-static/blob/master/create_dynamic.py#L44)
```python
wrf_file = 'wrfout_d03_2017-06-15_chch' # this is one example WRF output I used. The output file can be provided upon request.
```

3. Specify the start and end time stamp as well as the update frequency of boundary conditions [Line 53 in `create_dynamic.py`](https://github.com/dqbuhtig/WRF-PALM-no-static/blob/master/create_dynamic.py#L53)

```python
dt_start    = datetime(2017, 6, 16, 0,)   #  start time in YYYY/MM/DD/HH format
dt_end      = datetime(2017, 6, 17, 0,)   #  end time in YYYY/MM/DD/HH format
interval    = 2                           #  define time interval of WRF output to be read for the coupling
ts          = '2hour'                     #  specify the update frequency of boundary conditions which will 
                                          #  show in the dynamic input filename
                                          #  this works as a reference in case the update frequency calculation went wrong
```

**Remark:**   
- here users can define the start/end time and update frequency depending on their own requirments. For this example case, I have hourly WRF output, but I only want the boundary conditions update every 2 hours. The `interval` is hence set to 2. If you have 10-minute WRF output and want hourly update in PALM, then `interval = 6` (simple math :).  
- When the update frequency is not a divisor of the total run time, the script will round up the last step of time to avoid automatic termination of PALM. For example, if one has a 24-hour PALM simulation (`end_time = 86400`) with 5-hour update of boundary conditions, the time stamp in the dynamic file should be 0s, 18000s, 36000s, 54000s, 72000s, 86400s (instead of 90000s because 86400s may reach the final time step in WRF output). If the 86400s step is not in the dynamic file, PALM will terminate at 72000s step.


4. Specify the depth of 8 soil layers [Line 60 in `create_dynamic.py`](https://github.com/dqbuhtig/WRF-PALM-no-static/blob/master/create_dynamic.py#L60)

```python
dz_soil = np.array([0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86]) # this is the default setup in PALM
```

When the `create_dynamic.py` is successfully executed, the dynamic file will be stored in `dynamic_files`. 





