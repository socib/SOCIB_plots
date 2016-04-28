#!/usr/bin/env python

import os
import glob
import numpy as np
import netCDF4 as netcdf
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib  import colors
from matplotlib  import mpl
from matplotlib.path import Path
import matplotlib.patches as patches


#------------------------------------------------
# Plot the finite-element mesh and the topography
#------------------------------------------------

os.system('clear')
basemap_resolution = 'l'

meshdir='/home/ctroupin/DIVA/Global/mesh/L3/'
figdir='/home/ctroupin/DIVA/Global/figures/'
figtype='.jpg'
figname='global_mesh_L3bbb'

meshfile = 'mesh.dat'
meshtopofile = 'meshtopo.dat'

# Create figure directory if neccesary
if not(os.path.exists(figdir)):
    os.makedirs(figdir)

# Select colormaps
cmap=plt.cm.spectral_r

# Region of interest
# (could be read from GridInfo.dat file)
lonmin=-180.0
lonmax=180.0
latmin=-90.0
latmax=90.0

lonmap=0
latmap=20

#m = Basemap(resolution=basemap_resolution,projection='ortho',lat_0=latmap,lon_0=latmap)


#------------------------------------------------------------------------------

# Load mesh information
datamesh = np.loadtxt(meshdir+meshtopofile)
nnodes = int(datamesh[0])
ninterfaces = int(datamesh[1])
nelements = int(datamesh[2])
ntotal = int(datamesh.sum())

# Load mesh nodes
meshnodes = np.genfromtxt(meshdir+meshfile,skip_footer=nelements+ninterfaces)
meshnodes = np.fromstring(meshnodes)

# Load mesh elements
meshelements = np.genfromtxt(meshdir+meshfile,skip_header=nnodes+ninterfaces)
meshelements = np.fromstring(meshelements)
meshelements = np.int_(meshelements)

# Extract node coordinates
xnode=meshnodes[np.arange(1,nnodes*3,3)]
ynode=meshnodes[np.arange(2,nnodes*3,3)]

i=np.transpose(range(0,nnodes));
i1=meshelements[np.arange(0,nelements*6,6)]-1
i2=meshelements[np.arange(2,nelements*6,6)]-1
i3=meshelements[np.arange(4,nelements*6,6)]-1

# Make the plot
fig=plt.figure()
ax = fig.add_axes([0.1,0.2,0.8,0.7])



for j in range(0,nelements):
    verts = [(xnode[i1[j]],ynode[i1[j]]),\
         (xnode[i2[j]],ynode[i2[j]]),\
         (xnode[i3[j]],ynode[i3[j]]),\
         (xnode[i1[j]],ynode[i1[j]]) ]
    path = Path(verts)
    patch = patches.PathPatch(path, facecolor='none',edgecolor='grey',lw=0.5)
    ax.add_patch(patch)

ax.set_xlim(lonmin,lonmax)
ax.set_ylim(latmin,latmax)
#ax.drawcoastlines(ax=ax,color='black')

##m.drawparallels(np.arange(latmin,latmax,dlat), linewidth=0,
##                    labels=[1, 0, 0, 0], fontname='Times New Roman',fontsize=16)
##m.drawmeridians(np.arange(lonmin,lonmax,dlon), linewidth=0,
##                    labels=[0, 0, 1, 0], fontname='Times New Roman',fontsize=16)

#m.fillcontinents()
#m.bluemarble()
#m.etopo()
#m.drawcountries()
#m.drawrivers()
#m.nightshade(date)  

plt.savefig(figdir+figname+figtype, dpi=300, facecolor='w', edgecolor='w',
         transparent=False, bbox_inches='tight', pad_inches=0.1)
#plt.show()
plt.close()
