#!/usr/bin/env python

# plot_sst_avhrr_caibex_CTD_track.py
#
# Plot the SST 
#
#------------------------------------------------------------------------------

import os
import glob 
import numpy as np
import math
import netCDF4 as netcdf
import matplotlib.pyplot as plt
#import gmtColormap
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors

#-------------
# User options
#-------------

plotdaily=1

# Resolution for coastline 
basemap_resolution = 'h'

# File and directory names
CTDdir='/home/ctroupin/DataOceano/CAIBEX_campaign/data/CTD/'
CTDfile='CTD_coordinates_transect.dat'
seasoardir='/home/ctroupin/DataOceano/CAIBEX_campaign/data/Seasor/processedNoHead/Track1_6/'
seasoarfile='seasoar1_6.dat'
datadir='/data_local/Satellite/AVHRR/NAR/SSTonly/4Caibex/'
databasename='2009*.nc'
avgfile='SST_avg_Caibex.nc'     # Obtained using compute_sst_avg_avhrr_caibex.py
figdir='/home/ctroupin/Publis/201406_CAIBEX/figures/2015/'
coorddir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
coordfile='coord_NAR_SST.nc'

figtitleavg='2009-08-15 $-$ 2009-09-05'
fignameavg='SST_AVHRR_NAR_avg_colorbar'

# Figure extension
figtype1='.png'
figtype2='.eps'

# Colormap
##cdict= gmtColormap.gmtColormap('sst','/home/ctroupin/Software/Python/GMT_colormaps')
##cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
cmap = plt.cm.spectral

# Compute min and max values 
vmin=17.0
vmax=23.5
norm = colors.Normalize(vmin=vmin,vmax=vmax)
nlevels=150
levels2plot=np.arange(vmin,vmax+0.0001,(vmax-vmin)/(nlevels-1))
newticks=np.arange(np.floor(vmin),np.ceil(vmax)+0.0001,1.0)
newticks=newticks.round(2)
    
# Region of interest        
lonmin=-13.0
lonmax=-8.9
latmin=29.0
latmax=33.0

# Cape positions
lonSim,latSim=-(9.+(51./60.)),(31.+24./60.)
lonBed,latBed=-(9.+(16./60.)),(32.+33./60.)
lonGhir,latGhir=-(9.+(53.+20./60.)/60.),30.+(37.+49./60.)/60.
# text position
lonGhir2,latGhir2=-9.8,30.7
lonSim2,latSim2=-9.7,31.4
lonBed2,latBed2=-9.2,32.45

# Limits in the matrices (determined emperically)
imin,imax,jmin,jmax=450,800,1500,1800

# Spacing between x/y labels
dlon = 1
dlat = 1

# Unit of the variable to plot
unitname = ' '


m = Basemap(projection='merc',llcrnrlon=lonmin,llcrnrlat=latmin+0.1,\
                urcrnrlon=lonmax,urcrnrlat=latmax,  \
                lat_ts=0.5*(latmin+latmax),\
                resolution=basemap_resolution)

#------------------------------------------------------------------------------------


# Create figure directory if necessary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)

# Loop on files
filelist = sorted(glob.glob(datadir+databasename))
nfiles = len(filelist)

# Load CTD coordinates
coordCTD=np.loadtxt(CTDdir+CTDfile)
latCTD=coordCTD[:,0]+coordCTD[:,1]/60.+coordCTD[:,2]/3600.
lonCTD=-1*(coordCTD[:,3]+coordCTD[:,4]/60.+coordCTD[:,5]/3600.)
    
# Load coordinates from coordinate file
nc=netcdf.Dataset(coorddir+coordfile)
lon = nc.variables['lon'][imin:imax,jmin:jmax]
lat = nc.variables['lat'][imin:imax,jmin:jmax]
nc.close()

### Load SeaSoar data
lons,lats=np.loadtxt(seasoardir+seasoarfile,unpack=True,usecols=(0,1,))

# Mask bad values
lon[np.logical_or((lon<lonmin),(lon>lonmax))]=np.nan
lat[np.logical_or((lat<latmin),(lat>latmax))]=np.nan

xs,ys=m(lons,lats)


##print lon.min()
##print lon.max()
##print lat.min()
##print lat.max()

# For the capes and texts
xSim,ySim=m(lonSim,latSim)
xBed,yBed=m(lonBed,latBed)
xGhir,yGhir=m(lonGhir,latGhir)
# text position
xGhir2,yGhir2=m(lonGhir2,latGhir2)
xSim2,ySim2=m(lonSim2,latSim2)
xBed2,yBed2=m(lonBed2,latBed2)


meridians=np.arange(lonmin,lonmax+0.001,dlon)
parallels=np.arange(latmin,latmax+0.001,dlat)

# Projection
x,y = m(lon,lat)
xCTD,yCTD = m(lonCTD,latCTD)

NN=5
if plotdaily==1:
    for SSTfile in filelist:
        filename= os.path.basename(SSTfile)
        figname=filename.split(".")[0]
        figtitle=str(figname[0:4])+'-'+str(figname[4:6])+'-'+str(figname[6:8])
                                                        
        print figname

        # Load data from current file
        nc=netcdf.Dataset(SSTfile)
        field = nc.variables['sst'][0,imin:imax,jmin:jmax]
        nc.close()

        # Select sub-region
        field = field.squeeze()

        # Mask land values
        field = np.ma.masked_array(field,field==field.fill_value)
        field = np.ma.masked_array(field,field>=23.5)


        # Make the plot
        fig=plt.figure()
        ax = fig.add_subplot(111)
        m.ax=ax

        # Add CTD coordinates
        m.plot(xCTD,yCTD,'ko',ms=2,zorder=5)
        # Seasoar tracks
        m.plot(xs[::NN],ys[::NN],'k--',color='0.15',zorder=3)

        # Add grid, coastline and continent
        m.drawcoastlines(ax=ax,zorder=4)
        m.fillcontinents(color='.8', ax=ax,zorder=3)
        m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,zorder=2)
        m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

        #pcm=m.pcolormesh(x,y,field,cmap=cmap,norm=norm,zorder=2)
        pcm=m.scatter(x,y,s=2.5,c=field,cmap=cmap,norm=norm,edgecolor='None',zorder=2,marker='s')
        # Add the colorbar
        cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02,extend='both')
        cbar.set_ticks(newticks)

        xt,yt=m(-8.8,33.1)
        plt.text(xt,yt,r'($^{\circ}$C)',fontsize=18)
        plt.title(figtitle,fontsize=24)

        # Add capes position
            
        m.plot(xSim,ySim,'ko',ms=5,color='k',zorder=5)
        plt.text(xSim2,ySim2,'CS')
        m.plot(xGhir,yGhir,'ko',ms=5,color='k',zorder=5)
        plt.text(xGhir2,yGhir2,'CG')
        m.plot(xBed,yBed,'ko',ms=5,color='k',zorder=5)
        plt.text(xBed2,yBed2,'CB')  

        
    ##    ### Export figure and display it
        plt.savefig(figdir+figname+figtype1, dpi=300, facecolor='w', edgecolor='w',
                     transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.savefig(figdir+figname+figtype2, dpi=300, facecolor='w', edgecolor='w',
                     transparent=False, bbox_inches='tight', pad_inches=0.1)

        #plt.show()
        plt.close()

# ------------------
# Plot average SST
# ------------------

                                               
# Load data from current file
nc=netcdf.Dataset(datadir+avgfile)
field = nc.variables['SST'][:]
lon = nc.variables['longitude'][:]
lat = nc.variables['latitude'][:]
nc.close()

# Mask land values
field = np.ma.masked_array(field,field==field.fill_value)
field = np.ma.masked_array(field,field>=23.5)

# Make the plot
fig=plt.figure()
ax = fig.add_subplot(111)
m.ax=ax

# Add CTD coordinates
m.plot(xCTD,yCTD,'ko',ms=2,zorder=5)
# Seasoar tracks
m.plot(xs[::NN],ys[::NN],'k--',color='0.15',zorder=3)

# Add grid, coastline and continent
m.drawcoastlines(ax=ax,zorder=4)
m.fillcontinents(color='.8', ax=ax,zorder=3)
m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,zorder=2)
m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

pcm=m.pcolormesh(x,y,field,cmap=cmap,norm=norm,zorder=2)
#pcm=m.scatter(x,y,s=2.5,c=field,cmap=cmap,norm=norm,edgecolor='None',zorder=2,marker='s')
# Add the colorbar
cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02,extend='both')
cbar.set_ticks(newticks)

xt,yt=m(-8.8,33.1)
plt.text(xt,yt,r'($^{\circ}$C)',fontsize=18)
plt.title(figtitleavg,fontsize=24)

# Add capes position
    
m.plot(xSim,ySim,'ko',ms=5,color='k',zorder=5)
plt.text(xSim2,ySim2,'CS')
m.plot(xGhir,yGhir,'ko',ms=5,color='k',zorder=5)
plt.text(xGhir2,yGhir2,'CG')
m.plot(xBed,yBed,'ko',ms=5,color='k',zorder=5)
plt.text(xBed2,yBed2,'CB')  


##    ### Export figure and display it
plt.savefig(figdir+fignameavg+figtype1, dpi=300, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)
plt.savefig(figdir+fignameavg+figtype2, dpi=300, facecolor='w', edgecolor='w',
             transparent=False, bbox_inches='tight', pad_inches=0.1)

#plt.show()
plt.close()

