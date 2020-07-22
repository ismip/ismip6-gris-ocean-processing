# submarine_melt_implementation
Scripts relevant to the submarine melt implementation described in Slater et al. 2020, The Cryosphere (https://tc.copernicus.org/articles/14/985/2020/)

## EN4_ISMIP6_highres.mat
Estimated 'present day' ocean conditions, based on the EN4 dataset. Used to bias-correct the CMIP output in gen_projections_v2.m 

## gen_projections_v2.m
Takes CMIP ocean output and processes into sector by sector profiles, saves as ocean_extrap .mat files

## ocean_extrap .mat files
Contain sector ocean profiles for use as input to ComputeTF.m

## ComputeTF.m
Does extrapolation of ocean properties into fjords taking account of bathymetry, and saves output in netcdf format

## latlon2utm.m
Coordinate conversion from lat-lon (as provided in CMIP output) to UTM coordinates used in BedMachine
