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
import gmtColormap
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
from matplotlib.colors import LogNorm
from datetime import date

# Clean 
os.system('clear')

doplot=1

#-------------
# User options
#-------------

# Resolution for coastline 
basemap_resolution = 'h'

# File and directory names
CTDdir='/home/ctroupin/ULPGC/CAIBEX_campaign/data/CTD/'
CTDfile='CTD_coordinates_transect.dat'
datadir='/media/ctroupin/Iomega_HDD/DataOceano/Satellite/MODIS/data/daily/chlorophyll/'
databasename='A2009238_*.nc'
figdir='/home/ctroupin/Publis/CAIBEX_CSR2012/figures/2014/'
#figdir='/media/ctroupin/Iomega_HDD/DataOceano/Satellite/MODIS/figures/daily/chloro_python/'
coorddir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
coordfile='coord_NAR_SST.nc'

# Figure extension
figtype1='.png'
figtype2='.eps'

# Colormap
#cdict= gmtColormap.gmtColormap('sst','/home/ctroupin/Software/Python/GMT_colormaps')
#cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
#cmap=plt.cm.YlGnBu_r

cmap=plt.cm.gist_rainbow_r

# Compute min and max values 
vmin=0.029
vmax=3.1
norm = LogNorm(vmin=vmin,vmax=vmax)
nlevels=150
#newticks=np.arange(np.floor(vmin),np.ceil(vmax)+0.0001,1.0)
#newticks=newticks.round(2)

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
lonGhir2,latGhir2=-9.8,30.75
lonSim2,latSim2=-9.75,31.4
lonBed2,latBed2=-9.3,32.45

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
xt,yt=m(-9.500,33.1)

if doplot==1:
    
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
        m.plot(xCTD,yCTD,'ko',ms=2)

        

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
        
        # Export figure and display it
        plt.savefig(figdir+figname+'_V2', dpi=300, facecolor='w', edgecolor='w',
                     transparent=False, bbox_inches='tight', pad_inches=0.1)

        # Add capes position
            
        m.plot(xSim,ySim,'ko',ms=5,color='k'); 
        plt.text(xSim2,ySim2,'C.S.')
        m.plot(xGhir,yGhir,'ko',ms=5,color='k'); 
        plt.text(xGhir2,yGhir2,'C.G.')
        m.plot(xBed,yBed,'ko',ms=5,color='k'); 
        plt.text(xBed2,yBed2,'C.B.')

        plt.savefig(figdir+figname+'_capes_V2', dpi=300, facecolor='w', edgecolor='w',
                     transparent=False, bbox_inches='tight', pad_inches=0.1)
        

#        plt.show()
        plt.close()

