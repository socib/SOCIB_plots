#!/usr/bin/env python
# -------------------------------------------------------------------------

from matplotlib import colors
import datetime
import load_functions
import numpy as np
import matplotlib.pyplot as plt
import cf
from osgeo import gdal
from matplotlib._png import read_png
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
import os

doplotsst, doplotstream, plotsalinity = 0, 1, 1

# Regions of interest and box for the experiment
dlon, dlat = 1.0, 1.0
coordinates = np.array((0.5, 6.50001, 36.9, 41.))
coordinates2 = np.array((-2.00, 8.001, 35.5, 42.))
valex = -999.

year, month, day = 2016, 7, 26

yyyy, mm, dd = str(year), str(month).zfill(2), str(day).zfill(2)

figdate = yyyy + mm + dd
figtitle =  datetime.datetime(year, month, day).strftime("%d %B %Y")
figname = "Turtle_Altimetry_visible"+ figdate

romsfile = ("http://thredds.priv.socib.es/thredds/dodsC/operational_models/"
            "oceanographical/hydrodynamics/wmop/latest.nc")

#turtlefile = "http://thredds.priv.socib.es/thredds/dodsC/animal/turtle/turtle_mel-alk_ttrk014/L1/2015/08/dep0001_turtle-mel_alk-ttrk014_L1_2015-08-12.nc"

turtlefile1 = ("http://thredds.priv.socib.es/thredds/dodsC/animal/turtle/"
               "turtle_pixel-scb_splash001/L1/2016/06/"
               "dep0001_turtle-pixel_scb-splash001_L1_2016-06-22.nc")
turtlefile2 = ("http://thredds.priv.socib.es/thredds/dodsC/animal/turtle/"
               "turtle_eddy-scb_splash003/L1/2016/07/"
               "dep0001_turtle-eddy_scb-splash003_L1_2016-07-13.nc")

altimetryfile = ("http://thredds.priv.socib.es/thredds/dodsC/satellite/"
                 "altimetry/aviso/madt/L4/{0}/{1}/nrt_med_allsat_madt_uv_"
                 "{0}{1}{2}_{0}{1}{2}.nc.gz").format(yyyy, mm, dd)
altimetryfile2 = ("http://thredds.priv.socib.es/thredds/dodsC/satellite/"
                 "altimetry/aviso/madt/L4/{0}/{1}/nrt_med_allsat_madt_h_"
                 "{0}{1}{2}_{0}{1}{2}.nc.gz").format(yyyy, mm, dd)

turtlelogo = "/home/ctroupin/Presentations/figures4presentations/logo/turtle.png"
figdir = '/home/ctroupin/Pictures/SOCIB/'

visiblefile = "/data_local/Satellite/Visible/nasa-worldview-2016-07-27.tiff"

gtif = gdal.Open(visiblefile)
gtif.GetProjectionRef()

#Set the plot axis limits to the proper map coordinates.
arr = gtif.ReadAsArray()
trans = gtif.GetGeoTransform()
extent = (trans[0], trans[0] + gtif.RasterXSize*trans[1],
          trans[3] + gtif.RasterYSize*trans[5], trans[3])

vmin, vmax, dvar = 36.5, 38.5, 0.5

variable = 'sea_surface_salinity'
cmap = plt.cm.RdYlBu_r

norm = colors.Normalize(vmin=vmin, vmax=vmax)
normvel = colors.Normalize(vmin=0.0, vmax=0.5)
bounds = np.arange(vmin, vmax + .0001, dvar)

def add_logo_on_map(imagepath, ax, position, zoom, zorder):
    logo2plot = read_png(imagepath)
    imagebox = OffsetImage(logo2plot, zoom=zoom)
    # coordinates to position this image

    ab = AnnotationBbox(imagebox, position,
                        xybox=(0., 0.),
                        xycoords='data',
                        pad=0.0,
                        boxcoords="offset points")
    ab.zorder = zorder

    ax.add_artist(ab)

# Compute min and max values
sstmin, sstmax = 24., 28.
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .001, 1.0)

# Load turtle data
lonturtle1, latturtle1 = load_functions.load_turtle_coord(turtlefile1)
lonturtle2, latturtle2 = load_functions.load_turtle_coord(turtlefile2)

# Load model salinity
if plotsalinity:
    f = cf.read(romsfile)
    field2plot = f.select(variable)[0]
    lon = field2plot[0].coord('lon').array
    lat = field2plot[0].coord('lat').array
    field2plot = field2plot.array[0].squeeze()

    lon_ts, lat_ts = np.meshgrid(lon, lat)

# Load altimetry
lon_alti, lat_alti, u, v = load_functions.load_altimetry_aviso_uv(altimetryfile, coordinates2)
lon_alti2, lat_alti2, adt = load_functions.load_altimetry_aviso_adt(altimetryfile2, coordinates2)
adt = adt - adt.mean()

normalti = colors.Normalize(vmin=-0.15, vmax=0.15)


fig = plt.figure()
ax = plt.subplot(111)

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels

ax.set_xlim(coordinates[0], coordinates[1])
ax.set_ylim(coordinates[2], coordinates[3])

# Mask
# lonturtle = np.ma.masked_where(lonturtleQC >  2, lonturtle)
# latturtle = np.ma.masked_where(latturtleQC >  2, latturtle)

plt.imshow(arr[:3,:,:].transpose((1, 2, 0)), extent=extent, zorder=3,alpha=1)

plt.plot(lonturtle1[:], latturtle1[:], 'go-', zorder=5, markersize=1, linewidth=2.5)
add_logo_on_map(turtlelogo, ax, [lonturtle1[-2] , latturtle1[-2]], 1., zorder=6)

plt.plot(lonturtle2[:], latturtle2[:], 'go-', zorder=5, markersize=1, linewidth=2.5)
add_logo_on_map(turtlelogo, ax, [lonturtle2[-2] , latturtle2[-2]], 1., zorder=6)

# Add currents from altimetry
# ---------------------------

# Or streamfunction
if doplotstream:

    plt.streamplot(lon_alti, lat_alti, u, v, color=adt, norm=normalti, arrowsize=2,
               density=12, linewidth=.25, cmap=plt.cm.RdBu_r, zorder=5)


plt.title(figtitle, fontsize=20)
plt.savefig(os.path.join(figdir, figname), dpi=300, facecolor='w', edgecolor='w',
            transparent=False, bbox_inches='tight', pad_inches=0.1)

# plt.show()
plt.close()
