# process_CMIP
Convert and process netcdf CMIP output into .mat files

CMIP model output is downloaded as netcdf files. Usually this is one file per variable (potential temperature or salinity) per scenario (historical, RCP2.6 or RCP8.5). Unfortunately each download is slightly different and the format of the netcdf files once downloaded can vary a little. For example, some downloads included both the historical and RCP8.5 future scenario in the same file, some downloads involved one file per decade which then needed to be combined, some netcdf formatting had ocean depths negative and some had ocean depths positive, some had longitude in the range [-180,180] and some in the range [0,360] etc. Hence, the processing is individual to the model and scenario.

Aside from these complications, the processing reads the CMIP model output (lat, lon, depth, time, temperature and salinity), converts to thermal forcing and does the required spatial and depth-averaging, and saves as model and scenario-specific .mat files.

## ocean_projection.m
Matlab script to do the required processing - see in file for detailed workflow and comments - produces all the .mat files.

## .mat files
Contain thermal forcing time series for the indicated models and scenarios.

## latlon2utm.m
Coordinate conversion from lat-lon (as provided in CMIP output) to UTM coordinates used in BedMachine
