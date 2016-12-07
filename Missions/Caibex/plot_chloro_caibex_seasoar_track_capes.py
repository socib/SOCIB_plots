#!/usr/bin/env python
#
# plot_chloro_caibex_CTD_track_capes
#
# Plot the chlorophyll concentration, the CTD track and the cape names
#
# Check http://matplotlib.org/examples/pylab_examples/pcolor_log.html
# for log colorscales
#------------------------------------------------------------------------------

import os
import glob 
import numpy as np
import math
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
from matplotlib.colors import LogNorm
from datetime import date

# Clean 
os.system('clear')

doplotdaily=0
doplotmean=0
doplotL2=1

#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'h'

# File and directory names
CTDdir='/home/ctroupin/DataOceano/CAIBEX_campaign/data/CTD/'
CTDfile='CTD_coordinates_transect.dat'
seasoardir='/home/ctroupin/DataOceano/CAIBEX_campaign/data/Seasor/processedNoHead/Track1_6/'
seasoarfile='seasoar1_6.dat'
datadir='/data_local/Satellite/MODIS/data/daily/chlorophyll/'
databasename='T2009238_*.nc'
datadir2='/data_local/Satellite/MODIS/data/L2/'
datafile2='A2009238141000.nc'
datadiravg='/data_local/Satellite/MODIS/data/daily/avg/'
datafileavg='chloro_avg_Terra_Caibex.nc'
figdir='/home/ctroupin/Publis/201406_CAIBEX/figures/2015/'
coorddir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
coordfile='coord_NAR_SST.nc'

fignameavg='chloro_terra_avg'
figtitleavg='2009-08-15 $-$ 2009-09-05'


# Figure extension
figtype1='.png'
figtype2='.eps'

# Colormap
cmap=plt.cm.gist_rainbow_r

# Compute min and max values 
vmin=0.029
vmax=3.1
norm = LogNorm(vmin=vmin,vmax=vmax)
nlevels=150

newticks =np.array((0.003,0.01,0.03,0.1,0.3,1.0,3.0))
newlabels=np.array((0.03,0.1,0.3,1.0,3.0))

# Region of interest        
lonmin=-13.0
lonmax=-8.9
latmin=29.0
latmax=33.0

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

# Cape positions
lonSim,latSim=-(9.+(51./60.)),(31.+24./60.)
lonBed,latBed=-(9.+(16./60.)),(32.+33./60.)
lonGhir,latGhir=-(9.+(53.+20./60.)/60.),30.+(37.+49./60.)/60.
# text position
lonGhir2,latGhir2=-9.8,30.7
lonSim2,latSim2=-9.7,31.4
lonBed2,latBed2=-9.2,32.45

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
    
meridians=np.arange(lonmin,lonmax+0.001,dlon)
parallels=np.arange(latmin,latmax+0.001,dlat)

### Load SeaSoar data
lons,lats=np.loadtxt(seasoardir+seasoarfile,unpack=True,usecols=(0,1,))

### Mask bad values
##lon[np.logical_or((lon<lonmin),(lon>lonmax))]=np.nan
##lat[np.logical_or((lat<latmin),(lat>latmax))]=np.nan

xs,ys=m(lons,lats)


# For the capes and texts
xSim,ySim=m(lonSim,latSim)
xBed,yBed=m(lonBed,latBed)
xGhir,yGhir=m(lonGhir,latGhir)
# text position
xGhir2,yGhir2=m(lonGhir2,latGhir2)
xSim2,ySim2=m(lonSim2,latSim2)
xBed2,yBed2=m(lonBed2,latBed2)

# Load coordinates in the first file
file2load=filelist[0]
nc=netcdf.Dataset(file2load)
lon=nc.variables['lon'][:]
lat=nc.variables['lat'][:]
ttime=nc.variables['time'][:]
nc.close()

# Projections
llon,llat=np.meshgrid(lon,lat)
x,y = m(llon,llat)
xCTD,yCTD = m(lonCTD,latCTD)
xt,yt=m(-9.00,33.1)

NN=5

if doplotdaily==1:
    
    for chlorofile in filelist:
        filename= os.path.basename(chlorofile)
        figname=filename.split(".")[0]
        

        # Load data from current file
        nc=netcdf.Dataset(chlorofile)
        field = nc.variables['chloro'][:]
        ttime = ttime=nc.variables['time'][:]
        nc.close()

        

        # Extract date
        dd=date.fromordinal(ttime)
        month=dd.month
        day=dd.day
        figtitle=str(figname[1:5])+'-'+str(month)+'-'+str(day)
                                                        
        print figtitle
        
        # Select sub-region
        field = field.squeeze()
        field=np.flipud(field)
        
        # Mask land values
        field = np.ma.masked_array(field,np.isnan(field))

        # Make the plot
        fig=plt.figure()
        ax = fig.add_subplot(111)
        m.ax=ax
       
        pcm=m.pcolormesh(x,y,field,cmap=cmap,norm=norm)

        # Add CTD coordinates
        m.plot(xCTD,yCTD,'ko',ms=2,zorder=5)
        # Seasoar tracks
        m.plot(xs[::NN],ys[::NN],'k--',color='0.15',zorder=3)

        # Add grid, coastline and continent
        m.drawcoastlines(ax=ax)
        m.fillcontinents(color='.85', ax=ax,zorder=3)
        m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,zorder=2)
        m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

        # Add the colorbar
        cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02)
#        cbar.set_ticks(newticks)
        cbar.set_ticks(newticks)
        cbar.set_ticklabels(newlabels)
        
        plt.text(xt,yt,r'(mg m$^{-3}$)',fontsize=18)
        plt.title(figtitle,fontsize=24)
        
        # Add capes position
        m.plot(xSim,ySim,'ko',ms=5,color='k',zorder=5) 
        plt.text(xSim2,ySim2,'CS')
        m.plot(xGhir,yGhir,'ko',ms=5,color='k',zorder=5)
        plt.text(xGhir2,yGhir2,'CG')
        m.plot(xBed,yBed,'ko',ms=5,color='k',zorder=5)
        plt.text(xBed2,yBed2,'CB')

        plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',
                     transparent=False, bbox_inches='tight', pad_inches=0.1)
        

        #plt.show()
        plt.close()



# Load data from current file
nc=netcdf.Dataset(datadiravg+datafileavg)
field = nc.variables['chloro'][:]
lon  = nc.variables['longitude'][:]
lat  = nc.variables['latitude'][:]
nc.close()


# Mask land values
field = np.ma.masked_array(field,np.isnan(field))
field = np.flipud(field)


if doplotmean ==1:
    # Make the plot
    fig=plt.figure()
    ax = fig.add_subplot(111)
    m.ax=ax

    pcm=m.pcolormesh(x,y,field,cmap=cmap,norm=norm)

    # Add CTD coordinates
    m.plot(xCTD,yCTD,'ko',ms=2,zorder=5)
    # Seasoar tracks
    m.plot(xs[::NN],ys[::NN],'k--',color='0.15',zorder=3)

    # Add grid, coastline and continent
    m.drawcoastlines(ax=ax)
    m.fillcontinents(color='.85', ax=ax,zorder=3)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,zorder=2)
    m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

    # Add the colorbar
    cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02)
    #        cbar.set_ticks(newticks)
    cbar.set_ticks(newticks)
    cbar.set_ticklabels(newlabels)

    plt.text(xt,yt,r'(mg m$^{-3}$)',fontsize=18)
    plt.title(figtitleavg,fontsize=24)

    # Add capes position
    m.plot(xSim,ySim,'ko',ms=5,color='k',zorder=5) 
    plt.text(xSim2,ySim2,'CS')
    m.plot(xGhir,yGhir,'ko',ms=5,color='k',zorder=5)
    plt.text(xGhir2,yGhir2,'CG')
    m.plot(xBed,yBed,'ko',ms=5,color='k',zorder=5)
    plt.text(xBed2,yBed2,'CB')

    plt.savefig(figdir+fignameavg, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)


    #plt.show()
    plt.close()

if doplotL2 ==1:

    # Load L2 data
    with netcdf.Dataset(datadir2+datafile2, 'r') as nc:
      lon = nc.variables['Navigation_Data_longitude'][:]
      lat = nc.variables['Navigation_Data_latitude'][:]
      chla = nc.variables['Geophysical_Data_chlor_a'][:]
    # Mask
    NN=1
    chla[chla<0]=np.nan
    lonstart=200
    lonend=600
    latstart=0
    latend=1

    x,y=m(lon[latstart:-latend:NN,lonstart:-lonend:NN],lat[latstart:-latend:NN,lonstart:-lonend:NN])
    fig=plt.figure()
    ax = fig.add_subplot(111)
    m.ax=ax

    pcm=m.scatter(x,y,s=4,c=chla[latstart:-latend:NN,lonstart:-lonend:NN],
                  marker='s',cmap=cmap,norm=norm,edgecolor='none')

    # Add CTD coordinates
    m.plot(xCTD,yCTD,'ko',ms=2,zorder=5)
    # Seasoar tracks
    m.plot(xs[::NN],ys[::NN],'k--',color='0.15',zorder=3)

    # Add grid, coastline and continent
    m.drawcoastlines(ax=ax)
    m.fillcontinents(color='.85', ax=ax,zorder=3)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,zorder=2)
    m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

    # Add the colorbar
    cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02)
    #        cbar.set_ticks(newticks)
    cbar.set_ticks(newticks)
    cbar.set_ticklabels(newlabels)

    plt.text(xt,yt,r'(mg m$^{-3}$)',fontsize=18)
    plt.title('2009-8-27',fontsize=24)

    # Add capes position
    m.plot(xSim,ySim,'ko',ms=5,color='k',zorder=5) 
    plt.text(xSim2,ySim2,'CS')
    m.plot(xGhir,yGhir,'ko',ms=5,color='k',zorder=5)
    plt.text(xGhir2,yGhir2,'CG')
    m.plot(xBed,yBed,'ko',ms=5,color='k',zorder=5)
    plt.text(xBed2,yBed2,'CB')

    figname=datafile2[:-3]
    plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)


    #plt.show()
    plt.close()
