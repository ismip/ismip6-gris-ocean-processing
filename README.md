# ismip6-gris-ocean-processing
Matlab workflow for going from CMIP model output to projected retreat and submarine melt rates. A full description of the parameterisations, choices and workflow can be found in Slater et al. 2019, The Cryosphere (https://tc.copernicus.org/articles/13/2489/2019/) and Slater et al. 2020, The Cryosphere (https://tc.copernicus.org/articles/14/985/2020/)

# Workflow
### Convert and process netcdf CMIP output into .mat files
```process_CMIP/```

### Combine projected thermal forcing and subglacial runoff into projected submarine melt rate for each glacier in Greenland
```glaciers/```

### Contains the definition polygons for the 7 ice-ocean sectors
```final_region_def/```

### Takes projected submarine melt rates and converts into projected retreat
```final_projections/```
