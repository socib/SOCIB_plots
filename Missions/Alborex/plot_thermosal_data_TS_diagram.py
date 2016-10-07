#!/usr/bin/python
'''
Make a figure with a subplot showing the thermosalinograph data on one side
and the T-S diagram on the other subplot.
Requires the sea water toolbox in order to compute density values.
'''
import os
from matplotlib  import colors
from alborex_functions import *
import glob
from seabird.cnv import fCNV
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/ctroupin/Software/Python/seawater-3.3.2")
import seawater
import logging

def configure_logging():
    logger = logging.getLogger("alborex_TSdigram_logger")
    logger.setLevel(logging.DEBUG)
    # Format for our loglines
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # Setup console logging
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # Setup file logging as well
    fh = logging.FileHandler('/home/ctroupin/logs/alborex_TSdiagram_plot.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

logger = configure_logging()

# Domain configuration
dlon ,dlat = .25, .25
coordinates = np.array((-1.25, 0., 36.6, 37.75))
res = 'i'

# Files and directories
coastdir = "/home/ctroupin/IMEDEA/Cartex2014/data/coastline/"
coastfile = "coastline_cartex_f.txt"
thermosalfile = "http://thredds.socib.es/thredds/dodsC/research_vessel/thermosalinometer/socib_rv-scb_tsl001/L1/2014/05/dep0015_socib-rv_scb-tsl001_L1_2014-05-25.nc"
datadir = "/home/vessel/RTDATA/socib_rv/SCB-SBE9002/rawArchive/dep0007_socib-rv_scb-sbe9002_L1_2014-05-25/PROCESSED_SOCIB_half/"
figdir = "/home/ctroupin/public_html/Alborex/figures_20161006/"
figname = "alborex_thermosal_TSdiagram_norotate"

valex = 999

# for temperature
sstmin, sstmax = 18., 20.5
cmapsst = plt.cm.RdYlBu_r
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax+.1,1)

# for salinity
saltmin, saltmax = 36.5, 38.25
cmapsalt = plt.cm.RdYlBu_r
normsalt = colors.Normalize(vmin=saltmin, vmax=saltmax)
boundsalt = np.arange(saltmin, saltmax+.001, 0.5)

logger.info("Load coast")
loncoast, latcoast = alborex_load_coast(coastdir, coastfile, valex)

logger.info("Load thermosalinograph data")
lonThermo, latThermo, timeThermo, tempThermo, saltThermo, fluorThermo = alborex_load_thermosal(thermosalfile)

logger.info("Compute density for T-S couples")
density2plot = np.arange(26, 30,.25)
smin, smax, tmin, tmax = 36.5, 38.6, 12.5, 20.5
ds, dt = 0.05, 0.1
tvec = np.arange(tmin, tmax, dt)
svec = np.arange(smin, smax, ds)
ssvec, ttvec = np.meshgrid(svec, tvec)
density = seawater.eos80.dens0(ssvec, ttvec) - 1000.0

manual_locations = [(36.7, tmax-0.75),(36.9, tmax-0.5), (37.3, tmax-0.5), (37.6, tmax-0.5),
                    (37.9, tmax-0.5), (38.25, tmax-0.5), (smax-0.15, 19.5), (smax-0.15, 18.5), (smax-0.15, 17.5),
                    (smax-0.15, 16.5), (smax-0.15, 15.5), (smax-0.15, 14.5), (38.33, 12.7)]

m = Basemap(projection='merc', llcrnrlon=coordinates[0], llcrnrlat=coordinates[2],
                urcrnrlon=coordinates[1], urcrnrlat=coordinates[3],
                lat_ts=0.5 * (coordinates[2] + coordinates[3]), resolution=res)

fig = plt.figure(figsize=(16, 8))
ax = fig.add_subplot(1,2,2)
ax.set_anchor('N')
m.ax = ax

loncoast, latcoast = m(loncoast, latcoast)
lonThermo, latThermo = m(lonThermo, latThermo)

m.plot(loncoast,latcoast,'k-',lw=.5)

MM = 500
scat = m.scatter(lonThermo[:MM], latThermo[:MM], s=40, c=saltThermo[:MM], marker='o', edgecolor="None",
          zorder=2, cmap=cmapsalt, norm=normsalt)

cbar = fig.colorbar(scat, cmap=cmapsalt, orientation='horizontal', pad=0.025, aspect=15, shrink=1., extend='both')
cbar.set_ticks(boundsalt)

m.drawparallels(np.arange(coordinates[2], coordinates[3], dlat), linewidth=0.5,
                        labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
m.drawmeridians(np.arange(coordinates[0],coordinates[1], dlon), linewidth=0.5,
                        labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)
#m.drawmapscale(.5,35.5,1,37,50, barstyle='simple', units='km', fontsize=12,zorder=3)
os.path.join(figdir, figname)

m.fillcontinents(ax=ax, color='w', zorder=2)
#plt.title('Salinity Leg 1', fontsize=24)

ax = fig.add_subplot(1,2,1, adjustable='box', aspect=0.25)
ax.set_anchor('N')
filelist = sorted(glob.glob(datadir+'dalx1*cnv'))
for datafiles in filelist:
    profile = fCNV(datafiles)
    temperature = profile['TEMP']
    salinity = profile['PSAL']
    pressure = profile['DEPTH']
    potentialtemp = seawater.eos80.ptmp(salinity, temperature, pressure, pr=0)
    fluor = profile['flSP']
    plt.plot(salinity, potentialtemp, 'ko', ms=1)
#    plt.scatter(salinity, temperature, s=10, c=fluor, edgecolor='None')
cont = plt.contour(svec, tvec, density, levels=density2plot, colors='.65', linestyles=':')
plt.clabel(cont,inline=True, fmt='%1.2f', zorder=3, manual=manual_locations)

plt.title(' ', fontsize=24)
plt.ylabel('T ($^{\circ}$C)', fontsize=16)
plt.xlabel('S', fontsize=16)
plt.xlim(smin, smax)
plt.ylim(tmin, tmax)
#plt.title('Leg 1', fontsize=24)

plt.savefig(os.path.join(figdir, figname), dpi=300, facecolor='w', edgecolor='w',transparent=False, bbox_inches='tight', pad_inches=0.1)

# plt.show()
plt.close()
