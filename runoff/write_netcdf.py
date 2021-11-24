#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np
import time
import os

def write_netcdf(t, x, y, variables, filename, timeUnits='', description='', dims=None, overwrite=False):
   # Open the netcdf
   if os.path.isfile(filename):
      if not overwrite:
         print('WARNING: file {:s} exists. Use overwrite flag to replace. Exiting.'.format(filename))
         return
      else:
         os.remove(filename)

   dataset = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
   
   # Extra dimensions (in addition to time, x, y)
   if dims:
      for k,v in dims.items():
         dataset.createDimension(k,v)

   # Time
   dataset.createDimension('time', None)
   t_var = dataset.createVariable('time', np.float64, ('time',))
   t_var.units = timeUnits
   t_var[:] = t

   # Spatial coordinates
   dataset.createDimension('x', len(x))
   dataset.createDimension('y', len(y))
   x_var  = dataset.createVariable('x', np.float32, ('x',))
   x_var.units = 'm'; x_var.standard_name = 'projection_x_coordinate'; # TODO: make these input args
   y_var  = dataset.createVariable('y', np.float32, ('y',))
   y_var.units = 'm'; y_var.standard_name = 'projection_y_coordinate'; # TODO: make these input args
   x_var[:] = x
   y_var[:] = y

   # Variables
   for variable in variables:
      v_var = dataset.createVariable(variable['name'], np.float32, variable['dims'])
      v_var[:] = variable['z']
      v_var.units = variable['units']
      v_var.long_name = variable['long_name']
   
   # TODO: make these input args
   dataset.createDimension('nchar', 40)
   mapping_var = dataset.createVariable('mapping', 'S1', ('nchar',))
   mapping_var.grid_mapping_name = 'polar_stereographic'
   mapping_var.latitude_of_projection_origin = '90.'
   mapping_var.standard_parallel = '70.'
   mapping_var.straight_vertical_longitude_from_pole = '-45.'
   mapping_var.semi_major_axis = '6378137'
   mapping_var.inverse_flattening = '298.257223563'
   mapping_var.false_easting = '0.'
   mapping_var.false_northing = '0.'

   # General fields
   dataset.proj4 = '+init=epsg:3413'
   dataset.description = description
   dataset.history = 'Created ' + time.ctime(time.time())
   #dataset.source = 'netCDF4 python module tutorial'
   
   dataset.close()
   
