
# coding: utf-8

# **Goal:** plot a heatmap based on the platform locations.<br>
# We can use the tools developed in plot_TS_diagram_all.

# In[1]:

import plot_TS_diagram_all
import folium
import glob
import logging
import numpy as np
import os
import logging
from geopy.geocoders import Nominatim
from folium import plugins


# In[2]:

def configure_log():
    logdir = './log/'
    if not(os.path.exists(logdir)):
        os.mkdir(logdir)
    logfilename = os.path.join(logdir, 'Profiler_heatmap.log')

    # create logger 
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logfilename)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.CRITICAL)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info('+++ Start new run +++')
    return logger


# In[3]:

def main():
    
    datadir = "/data_local/DataOceano/CMEMS/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/profiler-glider"
    outputdir = "/home/ctroupin/public_html/LeafletMaps"

    logger = configure_log()
        
    if not(os.path.exists(outputdir)):
        os.mkdir(outputdir)
        logger.debug('Create output directory:')
        logger.debug(outputdir)

    datafilelist = sorted(glob.glob(os.path.join(datadir, '*.nc')))
    logger.info("Working on {0} files".format(len(datafilelist)))
    
    logger.info('Creating a new map')
    map_platform = folium.Map(location=[0.0, 0.0], zoom_start=2)

    coords_all = []
    # Loop on the files
    for datafiles in datafilelist:

        figbasename = os.path.basename(datafiles)[:-3]
        logger.info('Working on %s' %(figbasename))
        platformtype = figbasename.split('_')[-2]

        # Read the variables from the data file
        logger.info('Reading variables from netCDF file')
        time, time_units, lon, lat, depth, temp, psal = plot_TS_diagram_all.read_variables(datafiles)
                
	try:
            
 	    coords = [(lati, loni) for loni, lati in zip(lon, lat)]   
            # folium.PolyLine(coords, color="red", weight=3).add_to(map_platform)
            coords_all.append((lat.mean(), lon.mean()))
	    # folium.plugins.HeatMap(coords, radius=10).add_to(map_platform)
	except ValueError:
	    logger.error('Problem with coordinates')
    folium.plugins.HeatMap(coords).add_to(map_platform)
    logger.info('Saving leaflet map in %s directory' %(outputdir))
    map_platform.save(os.path.join(outputdir, 'test2.html'))
    logger.info('+++Finished+++')

# In[4]:

if __name__ == "__main__":
    main()

