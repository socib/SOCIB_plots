#!/usr/bin/env python
#
# plot_sst_avhrr_caibex_CTD_track.py
#
# Plot the SST averaged over period of interest 
#
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
from scipy import stats                     # to compute "nanmean", abstent in my numpy

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
datadir='/media/ctroupin/Iomega_HDD/DataOceano/Satellite/AVHRR/NAR/SSTonly/files4avg/'
databasename='2009*.nc'
figdir='/home/ctroupin/Publis/CAIBEX_CSR2012/figures/2014/'
coorddir='/home/ctroupin/DataOceano/Satellite/AVHRR/NAR/SSTonly/'
coordfile='coord_NAR_SST.nc'
figname='SST_AVHRR_NAR_avg'
figtitle='2009-08-15 $-$ 2009-09-05'

# Figure extension
figtype1='.png'
figtype2='.eps'

# Colormap
cdict= gmtColormap.gmtColormap('sst','/home/ctroupin/Software/Python/GMT_colormaps')
cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)


# Compute min and max values 
vmin=17.0
vmax=24.0
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

valex=-999.0
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

##print lon.min()
##print lon.max()
##print lat.min()
##print lat.max()


meridians=np.arange(lonmin,lonmax+0.001,dlon)
parallels=np.arange(latmin,latmax+0.001,dlat)

# Projections
x,y = m(lon,lat)
xCTD,yCTD = m(lonCTD,latCTD)

# For the capes and texts
xSim,ySim=m(lonSim,latSim)
xBed,yBed=m(lonBed,latBed)
xGhir,yGhir=m(lonGhir,latGhir)
# text position
xGhir2,yGhir2=m(lonGhir2,latGhir2)
xSim2,ySim2=m(lonSim2,latSim2)
xBed2,yBed2=m(lonBed2,latBed2)

# Allocate
SSTmatrix=np.zeros((nfiles,lon.shape[0],lon.shape[1]))

# Loop on the files to compute the average
tindex=0                   
for SSTfile in filelist:
    
    nc=netcdf.Dataset(SSTfile)
    field = nc.variables['sst'][0,imin:imax,jmin:jmax]
    nc.close()
    field[field<0.0]=np.nan
    SSTmatrix[tindex]=field
    tindex+=1
                   
SSTavg=stats.nanmean(SSTmatrix)


                                                    

if doplot==1:

    print figname

    # Mask land values
    #SSTavg = np.ma.masked_array(SSTavg,np.isnan(SSTavg))
    #SSTavg = np.ma.masked_array(SSTavg,SSTavg<10.)

    # Make the plot
    fig=plt.figure()
    ax = fig.add_subplot(111)
    m.ax=ax
   
    pcm=m.pcolormesh(x,y,SSTavg,cmap=cmap,norm=norm)

    # Add CTD coordinates
    #m.plot(xCTD,yCTD,'ko',ms=2)

    # Add grid, coastline and continent
    m.drawcoastlines(ax=ax)
    m.fillcontinents(color='.85', ax=ax,zorder=3)
    m.drawmeridians(meridians,labels=[0, 0, 0, 1],linewidth=0.5,fontsize=18,zorder=2)
    m.drawparallels(parallels,labels=[1, 0, 0, 0],linewidth=0.5,fontsize=18,zorder=2)

    xt,yt=m(-9.0,33.1)
    plt.title(figtitle,fontsize=24)
    
    # Export figure and display it
    plt.savefig(figdir+figname+figtype1, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)

    # Add capes position
        
    m.plot(xSim,ySim,'ko',ms=5,color='k'); 
    plt.text(xSim2,ySim2,'C.S.')
    m.plot(xGhir,yGhir,'ko',ms=5,color='k'); 
    plt.text(xGhir2,yGhir2,'C.G.')
    m.plot(xBed,yBed,'ko',ms=5,color='k'); 
    plt.text(xBed2,yBed2,'C.B.')  

    # Export figure and display it
    plt.savefig(figdir+figname+'_capes'+figtype1, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)

    # Add the colorbar and the units
    cbar=fig.colorbar(pcm,cmap=cmap,orientation='vertical',fraction=0.1,pad=0.02)
    cbar.set_ticks(newticks)
    plt.text(xt,yt,r'($^{\circ}$C)',fontsize=18)

    plt.savefig(figdir+figname+'_colorbar'+figtype1, dpi=300, facecolor='w', edgecolor='w',
                 transparent=False, bbox_inches='tight', pad_inches=0.1)
    
    #plt.show()
    plt.close()

