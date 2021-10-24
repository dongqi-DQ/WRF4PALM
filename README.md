# WRF4PALM v1.1 README

# Contents
1. What's new?
2. How to use WRF4PALM namelist?
3. Quick comparison between WRF and PALM dynamic driver

# what's new in v1.1?
- move from wrf-python to salem to modify RAM usage 
- use xarray instead of netCDF4 package to modify RAM usage
- apply multiprocessing to improve computation time 
- users now only need to edit namelist instead of editing the script  
- add surface variables (e.g. U10, V10, T2, and Q2) for surface NaN solver 
- read WRF projection info when locate PALM domain 
- allow users to specify the projection of PALM simulation
- geostrophic winds are estimated using geopotential height instead of pressure 

# How to use?
1. Download the entire source code to local
2. Provide your own WRF output
3. Edit namelist
4. Run WRF4PALM

# namelist 





# Remark  
water temperature in static driver
we might release a static driver generator using global data set from Google earth engine
E

# what's new?
- move from wrf-python to salem  
- use xarray instead of netCDF4 package to modify RAM usage
- apply multiprocessing to improve computation time 
- users now only need to edit namelist instead of editing the script  
- add surface variables (e.g. U10, V10, T2, and Q2) for surface NaN solver  

# How to use?
1. Download the entire source code to local
2. Provide your own WRF output
3. Edit namelist
4. Run WRF4PALM

# namelist 





# Remark  
water temperature in static driver
we might release a static driver generator using global data set from Google earth engine

