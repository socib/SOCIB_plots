#!/usr/bin/python

from netCDF4 import Dataset
import time, datetime
import shutil

def change_mode_format(mode):
    data_mode_dic = {"real time": 'R', "delayed time": 'D'}
    return data_mode_dic.get(mode, '')

def change_date_format(timein):
    timeout = datetime.datetime.strptime(timein, '%Y-%m-%dT%H:%M:%S+00:00').strftime('%Y-%m-%dT%I:%M:%SZ')
    return timeout

def seconds2juliandays(timein):
    date1950 = datetime.datetime(1950, 1, 1)
    date1970 = datetime.datetime(1970, 1, 1)
    deltatime = date1970.toordinal()-date1950.toordinal()
    timeout = timein / 86400. + deltatime
    return timeout

def change_file_format(filename):

  with Dataset(filename, 'r+', format='NETCDF4') as rootgrp:

    # Global attributes
    rootgrp.data_type = 'EGO glider time-series data'
    rootgrp.format_version = '1.0'
    rootgrp.platform_code = '99999'
    rootgrp.date_update = change_date_format(rootgrp.date_modified)
    rootgrp.data_mode = change_mode_format(rootgrp.data_mode)
    rootgrp.naming_authority = 'EGO'
    rootgrp.id = filename.split('.')[0]  # taken from file name... maybe something better to do
    rootgrp.source = "Glider observation"
    rootgrp.Conventions = "CF-1.4 EGO-1.0"
    rootgrp.geospatial_lat_min = str(rootgrp.geospatial_lat_min)
    rootgrp.geospatial_lat_max = str(rootgrp.geospatial_lat_max)
    rootgrp.geospatial_lon_min = str(rootgrp.geospatial_lon_min)
    rootgrp.geospatial_lon_max = str(rootgrp.geospatial_lon_max)
    rootgrp.time_coverage_start = change_date_format(rootgrp.time_coverage_start)
    rootgrp.time_coverage_end = change_date_format(rootgrp.time_coverage_end)

    # Rename TIME dimension
    rootgrp.renameDimension('time', 'TIME')

    # Rename Depth variable
    rootgrp.renameVariable('depth', 'DEPTH')

    # Rename TIME variable and add attributes
    rootgrp.renameVariable('time', 'TIME')
    TIME = rootgrp.variables['TIME']
    TIME.long_name = "Epoch time"
    TIME.units = "seconds since 1970-01-01T00:00:00Z"
    TIME.valid_min = "0."  # problem to be mentionned
    TIME.valid_max = "9000000000"  # problem to be mentionned
    TIME.QC_procedure, = "1"
    TIME.comment, = " "
    TIME.ancillary_variable = "TIME_QC"
    TIME.sdn_parameter_urn = "SDN:P01::ELTMEP01"
    TIME.sdn_uom_urn = "SDN:P061::UTBB"

    # Create LATITUDE variable with its attributes
    # (can not rename it because change of fill value).
    LATITUDE = rootgrp.createVariable('LATITUDE', 'f8', ('TIME',), fill_value=99999)
    LATITUDE.long_name = "Gps fixed latitude"
    LATITUDE.standard_name = "latitude"
    LATITUDE.units = 'degrees_north'
    LATITUDE.axis = "Y"
    LATITUDE.valid_min = "-90"
    LATITUDE.valid_max = "90"
    LATITUDE.QC_procedure = "1"
    LATITUDE.ancillary_variable = "POSITION_QC"
    LATITUDE.comment = " "
    LATITUDE.reference = "WGS84"
    LATITUDE.coordinate_reference_frame = "urn:ogc:crs:EPSG::4326"
    LATITUDE.sdn_parameter_urn = "SDN:P01::ALATZZ01"
    LATITUDE.sdn_uom_urn = "SDN:P061::DEGN"
    LATITUDE[:] = rootgrp.variables['latitude'][:]
    del rootgrp.variables['latitude']

    # Create LONGITUDE variable with its attributes
    # (can not rename it because change of fill value).
    LONGITUDE = rootgrp.createVariable('LONGITUDE', 'f8', ('TIME',), fill_value=99999)
    LONGITUDE.long_name = "Gps fixed longitude"
    LONGITUDE.standard_name = "longitude"
    LONGITUDE.units = 'degrees_east'
    LONGITUDE.axis = "X"
    LONGITUDE.valid_min = "-180"
    LONGITUDE.valid_max = "180"
    LONGITUDE.QC_procedure = "1"
    LONGITUDE.ancillary_variable = "POSITION_QC"
    LONGITUDE.comment = " "
    LONGITUDE.reference = "WGS84"
    LONGITUDE.coordinate_reference_frame = "urn:ogc:crs:EPSG::4326"
    LONGITUDE.sdn_parameter_urn = "SDN:P01::ALONZZ01"
    LONGITUDE.sdn_uom_urn = "SDN:P061::DEGE"
    LONGITUDE[:] = rootgrp.variables['longitude'][:]
    del rootgrp.variables['longitude']

    # Create JULD variable with its attributes
    # and to convert time to Julian day
    JULD = rootgrp.createVariable('JULD', 'f8', ('TIME',))
    JULD.long_name = 'Julian 1950 time'
    JULD.standard_name = 'time'
    JULD.units = 'days since 1950-01-01T00:00:00Z'
    JULD.valid_min = '0'
    JULD.valid_max = '90000'
    JULD.QC_procedure = 'test'
    JULD.comment = ' '
    JULD.axis = 'T'
    JULD.ancillary_variable = 'JULD_QC'
    JULD.sdn_parameter_urn = ' '
    JULD.sdn_uom_urn = 'SDN:P061::UTAA'
    JULD[:] = seconds2juliandays(rootgrp.variables['TIME'][:])

    # Create TIME_QC variable with its attributes
    TIME_QC = rootgrp.createVariable('TIME_QC', 'byte', ('TIME',), fill_value=-128)
    TIME_QC.long_name = 'Quality flag'
    TIME_QC.conventions = 'EGO reference table 2'
    TIME_QC.flag_values = '0, 1, 2, 3, 4, 5, 8, 9'
    TIME_QC.flag_meanings = 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed interpolated_value missing_value'
    TIME_QC[:] = 0

    # Create JULD_QC variable with its attributes
    JULD_QC = rootgrp.createVariable('JULD_QC', 'byte', ('TIME',), fill_value=-128)
    JULD_QC.long_name = 'Quality flag'
    JULD_QC.conventions = 'EGO reference table 2'
    JULD_QC.flag_values = '0, 1, 2, 3, 4, 5, 8, 9'
    JULD_QC.flag_meanings = 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed interpolated_value missing_value'
    JULD_QC[:] = 0

    # Create POSITION_QC variable with its attributes
    POSITION_QC = rootgrp.createVariable('POSITION_QC', 'byte', ('TIME',), fill_value=-128)
    POSITION_QC.long_name = 'Quality flag'
    POSITION_QC.conventions = 'EGO reference table 2'
    POSITION_QC.flag_values = '0, 1, 2, 3, 4, 5, 8, 9'
    POSITION_QC.flag_meanings = 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed interpolated_value missing_value'
    POSITION_QC[:] = 0


if __name__ == '__main__':
  import sys
  import os
  
  ifile = sys.argv[1]
  ofile = sys.argv[2]
  
  print ('Copying input file to new location:')
  print ('  input : ' + ifile)
  print ('  output: ' + ofile)
  odirs = os.path.dirname(ofile)
  if not os.path.exists(odirs):
    os.makedirs(odirs)
  shutil.copy2(ifile, ofile)
  
  print ('Reformatting new file...')
  change_file_format(ofile)
  
  print ('Done!')

