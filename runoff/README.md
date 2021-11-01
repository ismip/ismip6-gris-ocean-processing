# runoff
Scripts and output for processing of MAR runoff into basin-cumulative subglacial discharge time series and rasters.

## scripts
### runoff\_timeseries\_by\_basin.py
This python script reads MAR model output netcdfs, accumulates runoff per drainage basin at each time step, and outputs either a mat file with the time series of discharge or a netcdf with a time series of ice-sheet-wide rasters with pixel values representing discharge. The type of output can be toggled within the options in the script.

### runoff\_add\_glacier\_names.m
This matlab script reads the mat file output from runoff\_timeseries\_by\_basin.py and restructures the data with IDs and glacier names from Rignot and Mouginot (2012).

### Supporting scripts and files
**basins4highres\_xy.mat\_submergedarea.tif**: Raster containing the submerged frontal area of each terminus face. Necessary for converting units in runoff\_timeseries\_by\_basin.py

**tidewaterbasins\_rignotid.mat\_tidewaterbasins.tif**: Raster with basins delineated from outlet glacier locations identified in Rignot and Mouginot (2012).

**basins4highres\_xy.mat\_basins\_select.tif**: Raster with basins delineated automatically by finding all ice/ocean boundaries around the ice sheet.

**RignotMouginot2012\_TableS1\_GlacierList.csv**: Table S1 from Rignot and Mouginot (2012) with outlet glacier IDs, names, and locations.

**write_netcdf.py**: Script that writes netcdf with rasters of subglacial discharge.

**\*runoffBias.mat files**: Mat files with runoff biases per basin between MAR3.9 and RACMO2.3p2, estimated from 1995-2015). These are calculated in the runoff\_timeseries\_by\_basin.py script and saved for use because the processing takes a long time.


## output
### \*withGlacierIDs.mat files
Matlab structures (one per model and per scenario) containing projected subglacial runoff per tidewater glacier in Greenland. The data comes from running MAR3.9 regional climate model simulations over Greenland, which are themselves forced by the indicated CMIP model/scenario.


## References
Rignot, E., and J. Mouginot (2012), Ice flow in Greenland for the International Polar Year 2008–2009, Geophys. Res. Lett., 39, L11501, doi:10.1029/2012GL051634.