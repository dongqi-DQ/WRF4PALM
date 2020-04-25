# WRF-PALM coupler  
## for users who don't need to create new static files


This repository contains the following files:

- `create_cfg.py` - Python script to create a cfg file including PALM domain information  
- `create_dynamic.py` - the WRF-PALM coupling script  
- scripts in `util` - functions to use in the coupling  
- An example cfg file in `cfg_input`  
- An example dynamic input in `dynamic files`  

## Note:
- The scripts may only work on Python 3. 
- Vertical streching is currently not supported.  
- Users can downscale WRF to PALM when the PALM domain is smaller than one grid spacing of WRF, but this may not be very accurate.  
- Small shifts may apply when locating small domains using latitudes and logitudes.  
- When coupling to high resolution and large PALM domain, the coupling script may cost a lot of RAM. We hence recommend that .. 

