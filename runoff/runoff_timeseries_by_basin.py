#!/usr/bin/env python

import os, sys, argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
#from subprocess import call
import numpy as np
from netCDF4 import Dataset
#from pyproj import Proj, transform

import glob

#from scipy.interpolate import interp2d
from scipy import ndimage
from scipy.io import savemat, loadmat

import pickle

from osgeo import gdal, osr
from osgeo import gdalconst

import write_netcdf

import copy

import multiprocessing as mp
from progress.bar import Bar
import time

# ----------------
# Nested functions
# ----------------
def read_runoff(ncfilename,xVar,yVar,runoffVar,runoffCalculation):
#{{{
   # Read file
   ncfile = Dataset(ncfilename, 'r')

   x = ncfile.variables[xVar][:]
   y = ncfile.variables[yVar][:]
   lat = ncfile.variables['LAT'][:,:]
   lon = ncfile.variables['LON'][:,:]
   runoff = ncfile.variables[runoffVar]
   
   # Extract year
   year = int(os.path.basename(ncfilename).split('-')[-1].split('.')[0])

   if runoffCalculation == 'monthly':
      # Monthly means
      runoff_mean = runoff[:,:,:]
      for month in range(1,13):
         dt = (datetime.strptime('{:02d}/1/{:4.0f}'.format(month,year), '%m/%d/%Y') - timeEpoch).days
         t0 = datetime.strptime('{:02d}/01/{:4.0f}'.format(month,year), '%m/%d/%Y')
         t1 = t0 + relativedelta(months=1)
         time_bounds.append( np.array( [[(t0 - timeEpoch).days, \
                                         (t1 - timeEpoch).days]] ).T)

   if runoffCalculation == 'yearly-sum':
      # Annual mean
      runoff_mean = np.expand_dims(np.sum(runoff[:,:,:], axis=0), axis=0) # runoff / year
      dt = (datetime.strptime('07/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days
      time_bounds.append( np.array( [[(datetime.strptime('01/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days, \
                                      (datetime.strptime('12/31/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days]] ).T)

   if runoffCalculation == 'JJA-mean':
      # JJA mean
      runoff_mean = np.expand_dims(np.mean(runoff[5:8,:,:], axis=0), axis=0)
      dt = (datetime.strptime('07/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days
      time_bounds.append( np.array( [[(datetime.strptime('06/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days, \
                                      (datetime.strptime('08/31/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days]] ).T)

   return dt, time_bounds, x, y, runoff_mean, year
#}}}

def moving_average_downsample(x, y, window):
#{{{
   window = int(window)
   cumsum = np.cumsum(np.insert(y, 0, 0))
   y_mvavg = (cumsum[window:] - cumsum[:-window]) / float(window)
   y_mvavg = y_mvavg[::window]
   return y_mvavg
#}}}

def calculate_bias(basinNum, basinArrayBaseline, basinArray, runoffBaseline, runoff):
#{{{
   mask2d = basinArrayBaseline==basinNum
   mask3d = np.repeat(mask2d[np.newaxis,:,:], runoffBaseline.shape[0], axis=0)
   runoff_bl_basin = mask3d * np.flip(runoffBaseline, 1)
   mask2d = basinArray==basinNum
   mask3d = np.repeat(mask2d[np.newaxis,:,:], runoff.shape[0], axis=0)
   runoff_pd_basin = mask3d * np.flip(runoff, 1)
   runoff_bl_basin_timeseries = np.sum(runoff_bl_basin, axis=(1,2))
   runoff_pd_basin_timeseries = np.sum(runoff_pd_basin, axis=(1,2))

   # Yearly sums
   runoff_bl_basin_annual = 12. * moving_average_downsample(range(0,len(runoff_bl_basin_timeseries)), runoff_bl_basin_timeseries, 12)
   runoff_pd_basin_annual = 12. * moving_average_downsample(range(0,len(runoff_pd_basin_timeseries)), runoff_pd_basin_timeseries, 12)

   # Differences in annual means
   runoff_annual_diffs = runoff_pd_basin_annual - runoff_bl_basin_annual

   runoffBias = np.mean(runoff_annual_diffs)

   return runoffBias
#}}}

def runoff_unit_conversion(inputUnits, outputUnits):
#{{{
   # This assumes that the grid resolution is 1 km x 1 km
   if inputUnits == 'mmWE' and outputUnits == 'm3':
      unitConversion = 1000.
   elif inputUnits == 'mmWE yr-1' and outputUnits == 'kg s-1':
      unitConversion = 1000000./31557600.
   else:
      print('runoff unit conversion from ' + inputUnits + ' to ' + outputUnits + ' not supported!')
      return None
   return unitConversion
#}}}

def submerged_area(submergedAreaArray, basinArray, basinNum):
#{{{
   submergedAreaMasked = (basinArray==basinNum) * submergedAreaArray
   submergedArea = np.unique(submergedAreaMasked)
   submergedArea = submergedArea[submergedArea > 0.]
   if len(submergedArea) != 1:
      return None
   return submergedArea[0]
#}}}

def readRasterBandAsArray(filename, bandnum, rasterBandNoDataValue=None, clip_extent=None):
#{{{
   raster = gdal.Open(filename, gdalconst.GA_ReadOnly)
   rasterBand = raster.GetRasterBand(bandnum)
   rasterBandArray = rasterBand.ReadAsArray(0, 0, raster.RasterXSize, raster.RasterYSize).astype(np.float)

   if rasterBandNoDataValue is None:
      rasterBandNoDataValue = rasterBand.GetNoDataValue()
      if rasterBandNoDataValue is not None:
         rasterBandArray[rasterBandArray - rasterBandNoDataValue < np.finfo(float).eps] = np.nan
   else:
      rasterBandArray[rasterBandArray==rasterBandNoDataValue] = np.nan

   return rasterBandArray
#}}}

# -------------
# --- Setup ---
# -------------
netcdf_dirs = list()
# Select the climate model
MARver = 'MAR3.9'
CMIP = 'MIROC5' # 'MIROC5' | 'NorESM1' | 'CSIRO-Mk3.6' | 'HadGEM2-ES' | 'IPSL-CM5-MR' | 'ACCESS1.3' | 'CNRM-CM6' | 'UKESM1-CM6' | 'CNRM-ESM2'
rcp  = 'rcp85'  # 'histo' | 'rcp26' | 'rcp85' | 'ssp126' | 'ssp585'

approach = 'melt-rate' # 'retreat-rate' | 'melt-rate'

multiprocess_flag = False

if 'ssp' in rcp: ##{{{
   histo_end  = '2014'
   proj_start = '2015'
else:
   histo_end  = '2005'
   proj_start = '2006'
## }}}

# Selections for ISMIP6 retreat-rate (low-resolution) processing
if approach == 'retreat-rate': ##{{{
   inputUnits = 'mmWE' # 'mmWE yr-1' | 'mmWE'
   monthly_or_yearly = 'monthly' # 'yearly' | 'monthly'
   runoffCalculation = 'JJA-mean' # 'yearly-sum' | 'JJA-mean' | 'monthly'
   outputUnits = 'm3' # 'kg s-1' | 'm3'
   biasCorrection = False
   basinsRaster = 'tidewaterbasins_rignotid' # 'tidewaterbasins_rignotid' | 'basins4highres'
   basin_extrap_buffer_npixels = 0
   outputMatFile = True
   outputNetCDF  = False
##}}}

# Selections for ISMIP6 melt-rate (high-resolution) processing
if approach == 'melt-rate': ##{{{
   inputUnits = 'mmWE yr-1' # 'mmWE yr-1' | 'mmWE'
   monthly_or_yearly = 'monthly' # 'yearly' | 'monthly'
   runoffCalculation = 'yearly-sum' # 'yearly-sum' | 'JJA-mean' | 'monthly'
   outputUnits = 'kg s-1' # 'kg s-1' | 'm3'
   biasCorrection = True
   basinsRaster = 'basins4highres' # 'tidewaterbasins_rignotid' | 'basins4highres'
   basin_extrap_buffer_npixels = 3
   outputMatFile = False
   outputNetCDF  = True
##}}}

# MAR/RACMO directories
MARdir = 'MAR_DIR'
RACMOdir = 'RACMO_DIR'

# Output directory
outputDirectory = '.'

# --------------------------------------------------
# Should not need to change anything below this line
if rcp == 'histo':
   netcdf_dirs.append(MARdir + '/' + MARver + '/ISMIP6/GrIS/' + CMIP + '-' + rcp + '_1950_' + histo_end)
else:
   netcdf_dirs.append(MARdir + '/' + MARver + '/ISMIP6/GrIS/' + CMIP + '-histo_1950_' + histo_end)
   netcdf_dirs.append(MARdir + '/' + MARver + '/ISMIP6/GrIS/' + CMIP + '-' + rcp + '_' + proj_start + '_2100')
netcdf_filename_template = MARver.replace('MAR','MARv') + '-' + monthly_or_yearly + '-' + CMIP + '-*.nc'

runoffVar = 'RU'
runoffUnitConversion = runoff_unit_conversion(inputUnits, outputUnits)

tVar = None
xVar = 'x'
yVar = 'y'
xOverride = None
yOverride = None
projection = None # use native projection

timeUnits = 'days since 1900-1-1 00:00:00'
timeEpoch = datetime.strptime('1/1/1900', '%m/%d/%Y')
runoffLongName = 'cumulative basin runoff'

# Bias-correction for present-day runoff
netcdf_baseline = RACMOdir + '/runoff.1958-2017.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc'
xVar_baseline = 'lon'
yVar_baseline = 'lat'
runoffVar_baseline = 'runoffcorr'
runoff_baseline_time_indexes = range(444,684) # 01/15/1995 to 12/15/2014
netcdf_present_day_files = list()
for year in range(1995,2015):
   if year < int(proj_start):
      netcdf_present_day_files.append(MARdir + '/' + \
            MARver + '/ISMIP6/GrIS/' + CMIP + '-histo_1950_' + histo_end + '/' + MARver.replace('MAR','MARv') + '-' + monthly_or_yearly + '-' + CMIP + '-histo-{:4d}.nc'.format(year))
   else:
      netcdf_present_day_files.append(MARdir + '/' + \
            MARver + '/ISMIP6/GrIS/' + CMIP + '-' + rcp + '_' + proj_start + '_2100/' + MARver.replace('MAR','MARv') + '-' + monthly_or_yearly + '-' + CMIP + '-' + rcp  + '-{:4d}.nc'.format(year))

# Output file description
netcdf_description  = 'Cumulative runoff (' + runoffCalculation + \
      ') within Greenland marine-terminating outlet glacier drainage basins. Prepared for ISMIP6 by Denis Felikson (denis.felikson@nasa.gov).'

# Drainage basins
#{{{
if basinsRaster == 'tidewaterbasins_rignotid':
   clipfile = 'tidewaterbasins_rignotid.mat_tidewaterbasins.tif'
   submergedAreaRaster = None
if basinsRaster == 'basins4highres':
   clipfile = 'basins4highres_xy.mat_basins_select.tif'
   submergedAreaRaster = 'basins4highres_xy.mat_submergedarea.tif'
   outputUnits = outputUnits + ' m-2'

# All basins
basinfilter = 'basin=*'

# Parallel processing
nprocesses = 8
#}}}

# ------------------
# --- Processing ---
# ------------------
print("Setup:")
for netcdf_dir in netcdf_dirs:
      print(" input netcdf dir(s): " + netcdf_dir)
print(" clip file          : " + clipfile)
print(" output directory   : " + outputDirectory)
print(" basin filter       : " + basinfilter)
print(" ")

# Read NetCDF file
#{{{
print("Processing:")
print(" reading netcdf(s)")
dts = np.array([])
time_bounds = list()
runoffs = np.array([])

x = None
y = None
lat = None
lon = None

for netcdf_dir in netcdf_dirs:

   if multiprocess_flag:
      pool = mp.Pool(processes=nprocesses)
      results = [pool.apply_async(read_runoff, args=(ncfilename,xVar,yVar,runoffVar,runoffCalculation)) for ncfilename in sorted(glob.glob(netcdf_dir + '/' + netcdf_filename_template))]
      pool.close()
      pool.join()
      for p in results:
         dt, time_bounds, x, y, runoff_mean, year = p.get()
         if x is None:
            if xOverride is not None: x = xOverride
            if yOverride is not None: y = yOverride
         # Append
         dts = np.append(dts, dt)
         runoffs = np.append( runoffs, runoff_mean, axis=0 ) if runoffs.size else runoff_mean
   else:
      for ncfilename in sorted(glob.glob(netcdf_dir + '/' + netcdf_filename_template)):
         print(ncfilename)
         dt, time_bounds, x, y, runoff_mean, year = read_runoff(ncfilename,xVar,yVar,runoffVar,runoffCalculation)
         if x is None:
            if xOverride is None: x = ncfile.variables[xVar][:]
            else:                 x = xOverride
            if yOverride is None: y = ncfile.variables[yVar][:]
            else:                 y = yOverride
         # Append
         dts = np.append(dts, dt)
         runoffs = np.append( runoffs, runoff_mean, axis=0 ) if runoffs.size else runoff_mean

t = dts
runoff = runoffs
time_bounds = np.array(time_bounds)

# Geotransform
xStep = x[1]-x[0]
yStep = y[1]-y[0]
geoTransform = [x[0]-xStep/2, xStep, 0., y[-1]+yStep/2, 0., -yStep]
#}}}

runoff_nrows = int(runoff.shape[1])
runoff_ncols = int(runoff.shape[2])

# Create basin mask(s) from clip file
# Prelim checks {{{
if clipfile.split('.')[-1] == 'shp' and not checkForField(clipfile, basinfilter.split('=')[0]):
   print('ERROR in ' + __file__ + ': Attribute "' + basinfilter.split('=')[0] + '" not found in ' + clipfile)
   sys.exit()
#}}}

if clipfile.split('.')[-1] == 'shp': #{{{
   print(" ")
   print("\033[93mWARNING! This script uses the exterior of each polygon! \033[0m") 
   print(" ")

   if '*' in basinfilter:
      basinValues = getUniqueFieldValues(clipfile, basinfilter.split('=')[0])
   else:
      basinValues = [basinfilter.split('=')[1]]

   # Setup mask array(s)
   if projection:
      print("\033[93mERROR! This script is not set up to handle clipfiles that are shp but have a different projection than climate data! Exiting. \033[0m") 
      sys.exit()
   maskArray = np.zeros( (len(basinValues), runoff_nrows, runoff_ncols) )

   # Create basin mask
   print(" creating basin masks")
   driver = ogr.GetDriverByName("ESRI Shapefile")
   dataSource = driver.Open(clipfile, 0)
   layer = dataSource.GetLayer()

   for iValue, basinValue in enumerate(basinValues):
      AF = basinfilter.split('=')[0] + '=' + str(basinValue)
      layer.SetAttributeFilter(AF)
      for feature in layer:
         geometry = feature.GetGeometryRef()
         ring = geometry.GetGeometryRef(0)
         xClip = list()
         yClip = list()
         for ixy, xy in enumerate(ring.GetPoints()):
            xClip.append(xy[0])
            yClip.append(xy[1])
            
         maskArray_feature = clipImage(np.ones( (runoff_nrows, runoff_ncols) ), xClip, yClip, geoTransform)
         maskArray_feature[np.isnan(maskArray_feature)] = 0.
         maskArray[iValue,:,:] = maskArray[iValue,:,:] + maskArray_feature
#}}}
elif clipfile.split('.')[-1] == 'tif': #{{{
   # Resample mask to reslution of climate data
   print(" ")
   print("\033[93mWARNING! This script resamples the basin mask to the resolution of the climate data using nearest-neighbor! \033[0m") 
   print(" ")

   # Setup mask array(s)
   files = list()
   files.append(clipfile)
   if submergedAreaRaster: files.append(submergedAreaRaster)
   for f in files:
      f_resampled = os.path.basename(f).replace('.tif','_resampled.tif')
      if os.path.isfile(f_resampled):
         print("\033[93mWARNING! Overwriting " + f_resampled + "\033[0m") 
         os.remove(f_resampled)
      if projection:
         print("\033[93mWARNING! Reprojecting to specified projection of climate data:\n\t" + projection + "\033[0m") 
         cmdStr = 'gdalwarp {:s} {:s} -r near -t_srs "{:s}" -ot UInt32 -of GTiff -te {:f} {:f} {:f} {:f} -tr {:f} {:f}'.format( \
               f, f_resampled, \
               projection, \
               np.min(x)-xStep/2, np.min(y)-yStep/2, np.max(x)+xStep/2, np.max(y)+yStep/2, \
               geoTransform[1], -geoTransform[5])
      else: 
         cmdStr = 'gdal_translate {:s} {:s} -r near -ot UInt32 -of GTiff -projwin {:f} {:f} {:f} {:f} -tr {:f} {:f}'.format( \
               f, f_resampled, \
               np.min(x)-xStep/2, np.max(y)+yStep/2, np.max(x)+xStep/2, np.min(y)-yStep/2, \
               geoTransform[1], -geoTransform[5])
      os.system(cmdStr)

   basinArray         = readRasterBandAsArray(os.path.basename(clipfile).replace('.tif','_resampled.tif'), 1)
   if submergedAreaRaster: submergedAreaArray = readRasterBandAsArray(os.path.basename(submergedAreaRaster).replace('.tif','_resampled.tif'), 1)
   else: submergedAreaArray = np.ones(basinArray.shape)
   
   if '*' in basinfilter:
      basinValues = ['{:.0f}'.format(b) for b in sorted(set(basinArray[~np.isnan(basinArray)]))]
      basinValues.remove('0')
   else:
      basinValues = [basinfilter.split('=')[1]]

#}}}
else: #{{{
   print(" clipfile file type not supported")
   sys.exit()
# }}}

# Bias correction
if biasCorrection: ##{{{
   outputBiasFilename = '-'.join( (MARver, CMIP, rcp, basinsRaster, 'runoffBias') ) + '.mat'
   if os.path.isfile(outputBiasFilename):
      print('previously calculated biases found: ' + outputBiasFilename)
      runoffBiasDict = loadmat(outputBiasFilename)

   else:
      sys.stdout.write(' calculating bias correction\n')
      sys.stdout.flush()
      
      # Resample mask to reslution of climate data
      ncfile_bl = Dataset(netcdf_baseline, 'r')
      x_bl = ncfile_bl[xVar_baseline][:]; xStep_bl = x_bl[1]-x_bl[0]
      y_bl = ncfile_bl[yVar_baseline][:]; yStep_bl = y_bl[1]-y_bl[0]
      print(" ")
      print("\033[93mWARNING! This script resamples the basin mask to the resolution of the climate data using nearest-neighbor! \033[0m") 
      print(" ")
      cmdStr = 'gdal_translate {:s} {:s} -r near -ot UInt32 -of GTiff -projwin {:f} {:f} {:f} {:f} -tr {:f} {:f}'.format( \
         clipfile, os.path.basename(clipfile).replace('.tif','_resampled.tif'), \
         np.min(x_bl)-xStep_bl/2, np.max(y_bl)+yStep_bl/2, np.max(x_bl)+xStep_bl/2, np.min(y_bl)-yStep_bl/2, \
         xStep_bl, yStep_bl)
      os.system(cmdStr)
      basinArray_bl = readRasterBandAsArray(os.path.basename(clipfile).replace('.tif','_resampled.tif'), 1)
   
      runoff_pd = np.array([])
      # Read present-day "baseline" netcdf
      ncfile_bl = Dataset(netcdf_baseline, 'r')
      runoff_bl = ncfile_bl.variables[runoffVar_baseline][runoff_baseline_time_indexes,:,:]
      # Read present-day runoff over the same time period
      for netcdf_present_day_file in netcdf_present_day_files:
         ncfile_pd = Dataset(netcdf_present_day_file, 'r')
         runoff_pd = np.append(runoff_pd, ncfile_pd.variables[runoffVar][:,:,:].filled(0.), axis=0) if runoff_pd.size else ncfile_pd.variables[runoffVar][:,:,:].filled(0.)
   
      # Calculate biases per basin
      runoffBiasDict = dict() #np.nan * np.ones(len(basinValues))
      for iValue, basinValue in enumerate(basinValues):
         basinNum = float(basinValue)
         runoffBiasDict['basin'+basinValue] = calculate_bias(basinNum, basinArray_bl, basinArray, runoff_bl, runoff_pd)
   
      # Store biases
      savemat(outputBiasFilename, runoffBiasDict)
##}}}

# Extrapolate using nearest neighbor
basinSizes = {basinValue: np.sum(len(np.where(basinArray==float(basinValue))[0])) for basinValue in basinValues}
print(" buffering basins by {:2d} pixels".format(basin_extrap_buffer_npixels))
basinArrayExtrap = copy.deepcopy(basinArray)
for iteration in range(0,basin_extrap_buffer_npixels): #{{{
   bar = Bar('Iteration {:d}'.format(iteration), max=len(basinValues))
   for iValue, basinValue in enumerate(sorted(basinSizes, key=basinSizes.get, reverse=True)):
      basinNum = float(basinValue)
      # Buffer using the maximum filter
      basinArrayFilt = ndimage.filters.maximum_filter( np.where(basinArrayExtrap==basinNum, basinNum, 0.), size=3)
      # Remove pixels buffered into another basin
      basinArrayFilt[basinArrayExtrap > 0.] = 0.
      # iExtrapolate the basinArray
      basinArrayExtrap = np.where(basinArrayFilt==basinNum, basinNum, basinArrayExtrap)
      
      bar.next()
   bar.finish()
   
#}}}

# Initialize output dict
runoffDict = dict()
runoffDict['time'] = {'time': t[:], 'units': timeUnits}
runoffDict['runoff'] = {'units': outputUnits}
if biasCorrection:
   runoffDict['runoff']['bias'] = dict()
for basinValue in basinValues:
   runoffDict['runoff']['basin'+basinValue] = np.empty( len(t) )
   if biasCorrection:
      runoffDict['runoff']['bias']['basin'+basinValue] = np.empty( len(t) )

# Mask at each timestep for all basins {{{
print(" masking at each timestep for each basin")
runoffSums = np.nan * np.ones( (len(dts), len(basinValues)) )
#for itime in progressbar.progressbar( range(0, runoff.shape[0]) ):
for itime in range(0, runoff.shape[0]):
   print("  -> time: {:05d}".format(itime))
   runoff_timeSlice = runoff[itime,:,:]
   
   #for iValue, basinValue in progressbar.progressbar( enumerate(basinValues) ):
   for iValue, basinValue in enumerate(basinValues):
      basinNum = float(basinValue)
      #print("  -> basin {:s}".format(basinValue))
      runoffMasked = (basinArray==basinNum) * np.flipud(runoff_timeSlice)
      runoffSum = np.sum( runoffMasked )
   
      # Remove bias
      if biasCorrection:
         runoffSum = runoffSum - runoffBiasDict['basin'+basinValue]

      # Unit conversion
      runoffSum = runoffSum * runoffUnitConversion

      # Divide by submerged area
      if submergedAreaRaster:
         submergedArea = submerged_area(submergedAreaArray, basinArray, basinNum)
         if not submergedArea:
            import pdb; pdb.set_trace()
         runoffSum = runoffSum / submergedArea

      # If runoffSum < 0., set it to 0.
      if runoffSum < 0.:
         runoffSum = 0.

      # Store in dict
      runoffSums[itime,iValue] = runoffSum
      runoffDict['runoff']['basin'+basinValue][itime] = runoffSum
      
      if biasCorrection:
         runoffDict['runoff']['bias']['basin'+basinValue] = runoffBiasDict['basin'+basinValue]

print(' ')
# Output filename
outputFilename = '-'.join( (MARver, CMIP, rcp, runoffCalculation.replace('-','_'), basinsRaster) )

# Write mat file
if outputMatFile:
   print('writing mat file')
   if biasCorrection: outputFilename = outputFilename + '-biasCorrected'
   
   # Output mat file
   savemat(outputDirectory + '/' + outputFilename + '.mat', runoffDict)
#}}}

# For memory purposes, clear the runoff array
del runoff

# Write netcdf
if outputNetCDF:
   # Populate an array with the cumulative runoff
   print('writing nc file')
   runoffCumulativeBasins = np.nan * np.empty( (len(dts), basinArrayExtrap.shape[0], basinArrayExtrap.shape[1]) )
   for itime in range(0, len(dts)):
      for iValue, basinValue in enumerate(basinValues):
         basinNum   = float(basinValue)
         # Populate the runoff-by-basin array (cumulative within each basin)
         runoffCumulativeBasins[itime][np.where(basinArrayExtrap==basinNum)] = runoffSums[itime,iValue]
         
   # Write the netcdf
   dims = {'nv': 2}
   netcdf_variables = list()
   netcdf_variables.append({'name':'basin_runoff', 'z':np.flip(runoffCumulativeBasins, axis=1), 'dims':('time','y','x'), 'units':outputUnits, 'long_name':runoffLongName})
   netcdf_variables.append({'name':'time_bounds', 'z':time_bounds, 'dims':('time','nv'), 'units':timeUnits, 'long_name':'time bounds'})
   del runoffCumulativeBasins
   write_netcdf.write_netcdf(t, x, y, netcdf_variables, outputDirectory + '/' + outputFilename + '.nc', timeUnits=timeUnits, dims=dims, description=netcdf_description, overwrite=True)

