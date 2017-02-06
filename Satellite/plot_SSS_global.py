# Plot Sea Surface Salinity
#
# Data downloaded from http://podaac-opendap.jpl.nasa.gov/opendap/allData/aquarius/L3/mapped/V3/7day/SCID/2014/295/contents.html
# in NetCDF format (to avoid conversion)+
#

import glob
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from osgeo import gdal
import netCDF4 as netcdf
from matplotlib  import colors
import matplotlib.patheffects as path_effects
gdal.UseExceptions()

doplot=1
visibledir='/data_local/Satellite/Visible/'
#visiblefile='world.topo.200402.3x21600x10800.jpg'
visiblefile='world.topo.bathy.200402.3x5400x2700.jpg'
aquariusdir='/home/ctroupin/DataOceano/Satellite/Aquarius/data/'
aquariusfile='Q20142952014301.L3m_7D_SCID_V3.0_SSS_1deg.bz2.nc'
coastdir='/home/ctroupin/DataOceano/Coastlines/'
coastfile='medsea_coast1b.dat'
figdir='/home/ctroupin/DataOceano/Satellite/Aquarius/figures/'
figname='Q20142952014301_L3m_7D_SCID_V3_0_SSS_1deg'

cmap=plt.cm.RdYlBu_r
vmin,vmax,dvar=30,40,0.5
norm = colors.Normalize(vmin=vmin,vmax=vmax)
levels2plot=np.arange(vmin,vmax,dvar)
# region of interest
lonmin,lonmax,latmin,latmax,dlon,dlat=0,360,-80.,80.,0.2,0.2
lonmin2,lonmax2,latmin2,latmax2,dlon2,dlat2=-6.,35.,30.,48.,0.1,0.1

lon=np.arange(-180,180)
lat=np.arange(-89.5,89.5001,1.0)

# Load extracted coastline
valexc=-999.
lonc,latc=np.loadtxt(coastdir+coastfile,unpack=True)
loncoast=np.ma.masked_where(lonc==valexc,lonc)
latcoast=np.ma.masked_where(latc==valexc,latc)

with netcdf.Dataset(aquariusdir+aquariusfile) as nc:
    # Load coordinates
    SSS = nc.variables['_l3m_data'][:]
    SSS = np.flipud(SSS)
    
llon,llat=np.meshgrid(lon,lat)

if doplot==1:
    
    # create figure
    fig = plt.figure()
    ax = plt.subplot(111)

    # set up orthographic map projection with
    map = Basemap(projection='ortho',lat_0=5.,lon_0=-20,resolution='l')

    map.warpimage(image=visibledir+visiblefile,alpha=0.9,zorder=1)
    lonc,latc=map(llon,llat)
    lonc[lonc==1e+30]=np.nan
    latc[latc==1e+30]=np.nan

    #map.contourf(lonc,latc,SLA,levels2plot,cmap=cmap,norm=norm,zorder=2)
    map.pcolormesh(lonc,latc,SSS,cmap=cmap,norm=norm,zorder=2)

    
    plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',transparent=True, bbox_inches='tight', pad_inches=0.)
                
    #plt.show()
    plt.close

    fig=plt.figure()
    ax = fig.add_axes([0.1,0.2,0.7,0.5])
    
    map2 = Basemap(projection='merc',llcrnrlon=lonmin2,llcrnrlat=latmin2,\
            urcrnrlon=lonmax2,urcrnrlat=latmax2,  \
            lat_ts=0.5*(latmin2+latmax2),\
            resolution='h')

    lonc,latc=map2(llon,llat)
    loncoast,latcoast=map2(loncoast,latcoast)
    lonc[lonc==1e+30]=np.nan
    latc[latc==1e+30]=np.nan

    pcm=map2.pcolormesh(lonc,latc,SSS,cmap=cmap,norm=norm,zorder=2)
    map2.plot(loncoast,latcoast,'k-',lw=0.25)
    map2.fillcontinents()
    cbar=fig.colorbar(pcm,cmap=cmap,orientation='horizontal',fraction=0.1,pad=0.02)

    plt.savefig(figdir+figname+'_medsea', dpi=300, facecolor='w', edgecolor='w',transparent=True, bbox_inches='tight', pad_inches=0.)
                
    #plt.show()
    plt.close
