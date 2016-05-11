import numpy as np
from netCDF4 import Dataset
import time, datetime
import shutil

gliderdir = '/home/glider/data_dt/ideep02/20120509/netcdf/'
gliderfile = 'dep0003_ideep02_ime-sldeep002_L1_2012-05-09_data_dt.nc'
outputdir = '/home/ctroupin/SOCIB/Facilities/Glider/data/PerseusSOCIB_Glider/'
outputfilesuffix = '_EGOformat'

data_mode_dic = {"real time": 'R', "delayed time": 'D'}


def change_dataformat(timein):
    timeout = datetime.datetime.strptime(timein, '%Y-%m-%dT%H:%M:%S+00:00').strftime('%Y-%m-%dT%I:%M:%SZ')
    return timeout

def seconds2juliandays(timein):
    import datetime, time
    date1950 = datetime.datetime(1950, 1, 1)
    date1970 = datetime.datetime(1970, 1, 1)
    deltatime = date1970.toordinal()-date1950.toordinal()
    timeout = np.divide(timein,86400.)+deltatime
    return timeout

# Create new file name
outputfile = gliderfile[:-3] + outputfilesuffix + ".nc"
print outputfile

# Copy the glider file to another directory
shutil.copy2(gliderdir + gliderfile, outputdir + outputfile)

# Open netCDF file to modify names and attributes
rootgrp = Dataset(outputdir + outputfile, 'r+', format='NETCDF4')

# Dimension
rootgrp.renameDimension('time', 'TIME')


# Global attributes
rootgrp.data_type = 'EGO glider time-series data'
rootgrp.format_version = '1.0'
rootgrp.platform_code = '99999'
rootgrp.date_update = change_dataformat(rootgrp.date_modified)  # should be converted from rootgrp.date_modified
rootgrp.data_mode = data_mode_dic[rootgrp.data_mode]
rootgrp.naming_authority = 'EGO'
rootgrp.id = outputfile.split('.')[0]  # taken from file name... maybe something better to do
rootgrp.source = "Glider observation"
rootgrp.Conventions = "CF-1.4 EGO-1.0"
rootgrp.geospatial_lat_min = str(rootgrp.geospatial_lat_min)
rootgrp.geospatial_lat_max = str(rootgrp.geospatial_lat_max)
rootgrp.geospatial_lon_min = str(rootgrp.geospatial_lon_min)
rootgrp.geospatial_lon_max = str(rootgrp.geospatial_lon_max)

rootgrp.time_coverage_start = change_dataformat(rootgrp.time_coverage_start)
rootgrp.time_coverage_end = change_dataformat(rootgrp.time_coverage_end)

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


# Create LONGITUDE variable with its attributes
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


# Create JULD variable with its attributes   # Also need to convert time to Julian day
JULD = rootgrp.createVariable('JULD', 'f8', ('TIME',))
JULD.long_name = 'Julian 1950 time'
JULD.standard_name = 'time'
JULD.units = 'days since 1950-01-01T00:00:00Z'  # should be converted from existing value
JULD.valid_min = '0'
JULD.valid_max = '90000'
JULD.QC_procedure = 'test'
JULD.comment = ' '
JULD.axis = 'T'
JULD.ancillary_variable = 'JULD_QC'
JULD.sdn_parameter_urn = ' '
JULD.sdn_uom_urn = 'SDN:P061::UTAA'

# Create TIME_QC variable with its attributes

TIME_QC = rootgrp.createVariable('TIME_QC', 'byte', ('TIME',), fill_value=-128)
TIME_QC.long_name = 'Quality flag'
TIME_QC.conventions = 'EGO reference table 2'
TIME_QC.flag_values = '0, 1, 2, 3, 4, 5, 8, 9'
TIME_QC.flag_meanings = 'no_qc_performed good_data probably_good_data bad_data_that_are_potentially_correctable bad_data value_changed interpolated_value missing_value'

rootgrp.close()

print 'Finished processing file ' + str(gliderfile)