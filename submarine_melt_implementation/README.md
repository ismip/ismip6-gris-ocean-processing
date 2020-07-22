# submarine_melt_implementation
Scripts relevant to the submarine melt implementation described in Slater et al. 2020, The Cryosphere (https://tc.copernicus.org/articles/14/985/2020/)

## gen_projections_v2.m
Takes CMIP ocean output and processes into sector by sector profiles, saves as ocean_extrap .mat files

## ocean_extrap .mat files
Contain sector ocean profiles for use as input to ComputeTF.m

## ComputeTF.m
Does extrapolation of ocean properties into fjords taking account of bathymetry, and saves output in netcdf format
