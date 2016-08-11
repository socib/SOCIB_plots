#!/usr/bin/python
# coding: utf-8

# The goal is to compute and represent the Mixed Layer Depth (MLD) from the profiles acquired by the Argo profilers.<br>
# To do so, we will use the [python-oceans](https://github.com/pyoceans/python-oceans) packages, which has a function to compute the MLD.

# In[1]:

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import glob
import os
import seawater
import cmocean
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import FuncFormatter
import logging
import oceans
# %matplotlib inline
oceans.__version__


# ## Logging

# In[2]:

def configure_log():
    logdir = './log/'
    if not(os.path.exists(logdir)):
        os.mkdir(logdir)
    logfilename = os.path.join(logdir, 'computeMLD.log')
    # create logger 
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logfilename)
    fh.setLevel(logging.WARNING)
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
    return logger


# # Functions to read and plot

# In[3]:

def read_pressure_profiles(datafiles):
    '''
    Extract the temperature, salinity and coordinates from a netCDF file.
    Convert the depth to pressure if necessary.
    Keep only the good values (based on quality flags).
    '''
    with netCDF4.Dataset(datafiles) as nc:
        lon = nc.variables['LONGITUDE'][:]
        lat = nc.variables['LATITUDE'][:]
        try:
            pressure = nc.variables['PRES'][:]
        except KeyError:
            depth = nc.variables['DEPH'][:]
            pressure = seawater.eos80.pres(depth, np.tile(lat, (depth.shape[1], 1)).T)
            logger.info('Converting depth to pressure')
        try:
            temp = nc.variables['TEMP'][:]
            temp_qc = nc.variables['TEMP_QC'][:]
        except KeyError:
            logger.warning('Temperature not available in the profiles')
            temp = np.nan * np.ones_like(pressure)
            temp_qc = np.nan * np.ones_like(temp)
        try: 
            psal = nc.variables['PSAL'][:]
            psal_qc = nc.variables['PSAL_QC'][:]
        except KeyError:
            logger.warning('Salinity not available in the profiles')
            psal = np.nan * np.ones_like(temp)
            psal_qc = np.nan * np.ones_like(temp)
    temp, psal = apply_qc_variables(temp, temp_qc, psal, psal_qc)
    return lon, lat, abs(pressure), temp, psal


# In[4]:

def apply_qc_variables(temp, temp_qc, psal, psal_qc):
    '''
    Apply the quality flags (we keep only the good values)
    '''
    temp = np.ma.masked_where(temp_qc != 1, temp)
    psal = np.ma.masked_where(psal_qc != 1, psal)
    # chloro = np.ma.masked_where(chloro_qc > 3., chloro)
    return temp, psal


# In[5]:

def compute_mld_profile(pressure, temperature, salinity):
    if np.ma.isMaskedArray(pressure):
        pressure = pressure.compressed()
        MLD, idx_mld = oceans.ocfis.mld(salinity[:len(pressure)], temperature[:len(pressure)], 
                                    pressure, criterion='temperature')
    else:
        MLD, idx_mld = oceans.ocfis.mld(salinity, temperature, 
                                    pressure, criterion='density')
    return MLD


# In[6]:

def make_mld_plot(lon, lat, mld, figname, figtitle, **kwargs):
    m = Basemap(projection='robin', lon_0=0, resolution='c')
    lon, lat = m(lon, lat)
    fig = plt.figure(figsize=(10, 8))
    plt.title(figtitle)
    scat = plt.scatter(lon, lat, s=5, c=mld, edgecolor='None', alpha=0.75, **kwargs)
    cbar = plt.colorbar(extend='max', shrink=0.55)
    cbar.solids.set(alpha=1)
    cbar.set_label('MLD\n(m)', rotation=0, ha='left')
    cbar.set_ticks(range(0, 151, 25))
    m.fillcontinents(color='grey', zorder=3)
    # m.drawcoastlines(linewidth=0.25, zorder=4)
    # plt.show()
    plt.savefig(figname, dpi=300)
    plt.close()


# In[7]:

def make_mld_hist(mld, figname):
    mld = np.array(mld)
    mld = np.ma.masked_where(np.isnan(mld), mld)
    
    bins = np.linspace(0., 150., 16.)
    # bins = np.append(bins, (500, 750))
    plt.figure(figsize=(10,8))
    ax = plt.gca()
    plt.hist(mld, range=(0, 150.), bins=bins, 
             histtype='stepfilled', orientation="horizontal", color='k', 
             )
    plt.xlabel('Profile frequency')
    plt.ylabel('Mixed layer\ndepth (m)', rotation=0, ha='right')
    ax.invert_yaxis()
    # plt.show()
    plt.savefig(figname, dpi=300)
    plt.close()


# # Main loop

# In[10]:

def main():

    # logger = configure_log()
    period = '201605'
    figdir = "/home/ctroupin/Projects2/201501_InsTAC/Graphical_Material/2016_Q2/MLD"
    datadir = "/data_local/DataOceano/CMEMS/INSITU_GLO_NRT_OBSERVATIONS_013_030/monthly/profiler-glider/" + period
    
    if not(os.path.exists(figdir)):
        os.mkdir(figdir)
        logger.debug('Create figure directory')

    datafilelist = sorted(glob.glob(os.path.join(datadir, '*.nc')))
    nfiles = len(datafilelist)
    logger.info("Working on {0} files".format(nfiles))
    
    logger.info("Initialize lon, lat and mld list")
    lon_all = []
    lat_all = []
    mld_all = []
    
    # Loop on the files
    fcounter = 0
    for datafiles in datafilelist:

        fcounter += 1
        figbasename = os.path.basename(datafiles)[:-3]
        logger.info('Working on file {0} ({1}/{2})'.format(figbasename, fcounter, nfiles))
        platformtype = figbasename.split('_')[-2]
        
        # Check the platform type (PF = profiler or GL = glider)
        if platformtype == 'PF':
            logger.info('Reading variables from netCDF file')
            lon, lat, pressure, temp, psal = read_pressure_profiles(datafiles)
            nprofiles, ndepth = temp.shape
            
            # Check is temp and/or salt exist
            if (np.isnan(temp).all()) or (np.isnan(psal).all()):
                logger.info('No salinity or/and temperature to plot')
            else:
                # Check is the profiles have more than one depth
                if (ndepth > 3):
                    logger.debug('Computing MLD')
                    logger.debug('Loop on the %s profiles' %(nprofiles))
                
                    for i in range(0, nprofiles):
                        try:
                            mld = compute_mld_profile(pressure[i, :], temp[i, :], psal[i, :])
                            lon_all.append(lon[i])
                            lat_all.append(lat[i])
                            mld_all.append(mld)
                        except IndexError:
                            logger.error('arrays used as indices must be of integer (or boolean) type')
                else:
                    logger.info('Not enough depth levels to compute MLD')
                
                    
    plt.style.use('../../../stylefiles/socib.mplstyle')
    logger.info('Creating scatter plot')                
    make_mld_plot(lon_all, lat_all, mld_all, os.path.join(figdir, 'MLD_scatter_' + period), 'May 2016',
                  vmin=5, vmax=150., cmap=cmocean.cm.density)
    
    logger.info('Histogram') 
    make_mld_hist(mld_all, os.path.join(figdir, 'MLD_histogram' + period))


# In[11]:

if __name__ == '__main__':
    logger = configure_log()
    main()

