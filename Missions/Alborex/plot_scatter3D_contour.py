#!/usr/bin/python
#
# -------------------------------------------------------------------------

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
from alborex_functions import *

plottemp, plotsalinity, plotchloro = 1, 0, 0

# Plot configuration
dlon, dlat = 1., 1.
coordinates = np.array((-0.85, -0.25, 36.71, 37.19))
depthmin = 0.0
depthmax = 300.0
res = 'i'

# Files and directories
coastdir = '/home/ctroupin/IMEDEA/Cartex2014/data/coastline/'
coastfile = 'coastline_cartex_f.txt'
ctdfile = ('http://thredds.priv.socib.es/thredds/dodsC/research_vessel/ctd/',
           'socib_rv-scb_sbe9002/L1/2014/dep0007_socib-rv_scb-sbe9002_L1_2014-05-25.nc')
valex = 999
cmap = plt.cm.YlGnBu

tempCTDmin, tempCTDmax = 13., 19.
boundsCTDtemp = np.arange(tempCTDmin, tempCTDmax, 1.)
normCTDtemp = colors.Normalize(vmin=tempCTDmin, vmax=tempCTDmax)
cmapCTDtemp = plt.cm.spectral

loncoast, latcoast = alborex_load_coast(coastdir, coastfile, valex)

lonCTD, latCTD, depthCTD, tempCTD, chloroCTD = alborex_load_ctd(ctdfile)


depths = np.array((-600, 0, 100))
angles = np.array((55, -40))
props = dict(boxstyle='round', facecolor='white', alpha=0.99)


    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')

    logger.info("Plot the CTD profiles")
    ndepth = depthCTD.shape[1]
    for p in range(0, 20):
        scat3D = ax.scatter(lonCTD[p] * np.ones(ndepth), latCTD[p] * np.ones(ndepth),
                            -depthCTD[p, :], s=10, c=tempCTD[p, :], marker='o',
                            edgecolor='none', norm=normCTDtemp,
                            cmap=cmapCTDtemp, zorder=4)

    # Plot the glider profiles
    scat3Dglider1 = ax.scatter(longlider1, latglider1, -depthglider1, s=20,
                               c=temperatureglider1, edgecolor='none',
                               norm=normCTDtemp, marker='o',
                               cmap=cmapCTDtemp, zorder=4)
    scat3Dglider2 = ax.scatter(longlider2, latglider2, -depthglider2, s=20,
                               c=temperatureglider2, edgecolor='none',
                               norm=normCTDtemp, cmap=cmapCTDtemp, zorder=4)

    # Plot the surface CTD and glider tracks
    # ax.plot3D(lonCTD,latCTD,0,'--',color="0.45")
    # ax.plot3D(longlider2,latglider2,0,'--',color="0.25")
    # ax.plot3D(longlider1,latglider1,0,'--',color="0.25")

    # Add labels for the CTD
    # for p in range(0,40):
    #    ax.text(lonCTD[p],latCTD[p],15.,str(p),ha='center',zorder=6,bbox=props)

    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmap, orientation='vertical', pad=0.05,
                        aspect=15, shrink=0.8, norm=normCTDtemp, extend='both')
    cbar.set_label('$^{\circ}$C', fontname='Times New Roman',
                   rotation=0, ha='left', fontsize=18)
    cbar.set_clim(tempCTDmin, tempCTDmax)
    cbar.set_ticks(boundsCTDtemp)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(os.path.join(figdir, figname1 + '_gliderCTD'), dpi=300,
                facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close()

    logger.info("--- Coastal glider")

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')
    # Plot the glider profiles
    scat3Dglider2 = ax.scatter(longlider2, latglider2, -depthglider2,
                               s=20, c=temperatureglider2, edgecolor='none',
                               norm=normCTDtemp, cmap=cmapCTDtemp, zorder=4)

    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmap, orientation='vertical', pad=0.05,
                        aspect=15, shrink=0.8, norm=normCTDtemp, extend='both')
    cbar.set_label('$^{\circ}$C', fontname='Times New Roman', fontsize=18)
    cbar.set_clim(tempCTDmin, tempCTDmax)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname1 + '_glider1', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()

    logger.info("--- Deep glider")

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')

    # Plot the glider profiles
    scat3Dglider1 = ax.scatter(longlider1, latglider1, -depthglider1, s=5, c=temperatureglider1, edgecolor='none',
                               norm=normCTDtemp, cmap=cmapCTDtemp, zorder=4)

    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmap, orientation='vertical', pad=0.05, aspect=15, shrink=0.8,
                        norm=normCTDtemp,
                        extend='both')
    cbar.set_label('$^{\circ}$C', fontname='Times New Roman', fontsize=18)
    cbar.set_clim(tempCTDmin, tempCTDmax)
    cbar.set_ticks(boundsCTDtemp)
    cbar.ax.set_xticklabels(boundsCTDtemp, fontname='Times New Roman', fontsize=16)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname1 + '_glider2', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()

    # -----------------

if plotsalinity:

    depthCTD = np.ma.masked_outside(depthCTD, depthmax, abs(depthmin))
    depthglider1 = np.ma.masked_outside(depthglider1, depthmax, abs(depthmin))
    depthglider2 = np.ma.masked_outside(depthglider2, depthmax, abs(depthmin))

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')

    # Plot the CTD profiles
    ndepth = depthCTD.shape[1]
    for p in range(0, len(lonCTD)):
        scat3D = ax.scatter(lonCTD[p] * np.ones(ndepth), latCTD[p] * np.ones(ndepth), -depthCTD[p, :], s=15,
                            c=psalCTD[p, :], edgecolor='None', norm=normchloro, alpha=1, cmap=cmapCTDpsal, zorder=4)

    ##    # Plot the glider profiles
    scat3Dglider1 = ax.scatter(longlider1, latglider1, -depthglider1, s=5, c=chloroglider1, edgecolor='None',
                               norm=normCTDpsal, cmap=cmapCTDpsal, zorder=4)
    scat3Dglider2 = ax.scatter(longlider2, latglider2, -depthglider2, s=5, c=chloroglider2, edgecolor='None',
                               norm=normCTDpsal, cmap=cmapCTDpsal, zorder=4)
    #
    #    # Plot the surface CTD and glider tracks
    #    ax.plot3D(lonCTD,latCTD,0,'--',color="0.45")
    #    ax.plot3D(longlider2,latglider2,0,'--',color="0.25")
    #    ax.plot3D(longlider1,latglider1,0,'--',color="0.25")
    #
    #    ### Add labels for the CTD
    #    ##for p in range(0,40):
    #    ##    ax.text(lonCTD[p],latCTD[p],15.,str(p),ha='center',zorder=6,bbox=props)
    #
    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmapCTDpsal, orientation='vertical', pad=0.05, aspect=15, shrink=0.8, norm=normCTDpsal,
                        extend='both')
    cbar.set_clim(psalCTDmin, psalCTDmax)
    boundsCTDchloro = np.array((0.01, 0.03, 0.1, 0.3, 1.0))
    cbar.set_ticks(boundsCTDpsal)
    cbar.ax.set_yticklabels(boundsCTDpsal, fontname='Times New Roman', fontsize=16)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname2 + '_gliderCTD', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()

    # ----------------------------------------------

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')


    #    # Plot the glider profiles
    scat3Dglider1 = ax.scatter(longlider1, latglider1, -depthglider1, s=5, c=salinityglider1, edgecolor='None',
                               norm=normchloro, cmap=cmapCTDchloro, zorder=4)

    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmap, orientation='vertical', pad=0.05, aspect=15, shrink=0.8,
                        norm=normCTDpsal,
                        extend='both')
    cbar.set_clim(psalCTDmin, psalCTDmax)
    boundsCTDchloro = np.array((0.01, 0.03, 0.1, 0.3, 1.0))
    cbar.set_ticks(boundsCTDpsal)
    cbar.ax.set_yticklabels(boundsCTDpsal, fontname='Times New Roman', fontsize=16)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)hloro
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname2 + '_glider1', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()

    # ----------------------------------

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')

    # Plot the glider profiles
    scat3Dglider2 = ax.scatter(longlider2, latglider2, -depthglider2, s=5, c=salinityglider2, edgecolor='None',
                               norm=normCTDpsal, cmap=cmapCTDpsal, zorder=4)

    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmapCTDpsal, orientation='vertical', pad=0.05, aspect=15, shrink=0.8,
                        norm=normCTDpsal, extend='both')
    cbar.set_clim(psalCTDmin, psalCTDmax)
    cbar.set_ticks(boundsCTDpsal)
    cbar.ax.set_yticklabels(boundsCTDpsal, fontname='Times New Roman', fontsize=16)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname2 + '_glider2', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()

if plotchloro == 1:

    depthmin, depthmax, deltadepth = -300, 0, 50
    depths = np.array((depthmin, depthmax, deltadepth))

    depthCTD = np.ma.masked_outside(depthCTD, depthmax, abs(depthmin))
    depthglider1 = np.ma.masked_outside(depthglider1, depthmax, abs(depthmin))
    depthglider2 = np.ma.masked_outside(depthglider2, depthmax, abs(depthmin))

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')

    # Plot the CTD profiles
    ndepth = depthCTD.shape[1]
    for p in range(0, len(lonCTD)):
        scat3D = ax.scatter(lonCTD[p] * np.ones(ndepth), latCTD[p] * np.ones(ndepth), -depthCTD[p, :], s=15,
                            c=chloroCTD[p, :],
                            edgecolor='None', norm=normchloro, alpha=1, cmap=cmapCTDchloro, zorder=4)

    ##    # Plot the glider profiles
    scat3Dglider1 = ax.scatter(longlider1, latglider1, -depthglider1, s=5, c=chloroglider1, edgecolor='None',
                               norm=normchloro, cmap=cmapCTDchloro, zorder=4)
    scat3Dglider2 = ax.scatter(longlider2, latglider2, -depthglider2, s=5, c=chloroglider2, edgecolor='None',
                               norm=normchloro, cmap=cmapCTDchloro, zorder=4)
    #
    #    # Plot the surface CTD and glider tracks
    #    ax.plot3D(lonCTD,latCTD,0,'--',color="0.45")
    #    ax.plot3D(longlider2,latglider2,0,'--',color="0.25")
    #    ax.plot3D(longlider1,latglider1,0,'--',color="0.25")
    #
    #    ### Add labels for the CTD
    #    ##for p in range(0,40):
    #    ##    ax.text(lonCTD[p],latCTD[p],15.,str(p),ha='center',zorder=6,bbox=props)
    #
    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmap, orientation='vertical', pad=0.05, aspect=15, shrink=0.8, norm=normCTDchloro,
                        extend='both')
    cbar.set_label('$mg\,m^{-3}$', fontname='Times New Roman', fontsize=18)
    cbar.set_clim(chloroCTDmin, chloroCTDmax)
    boundsCTDchloro = np.array((0.01, 0.03, 0.1, 0.3, 1.0))
    cbar.set_ticks(boundsCTDchloro)
    cbar.ax.set_yticklabels(boundsCTDchloro, fontname='Times New Roman', fontsize=16)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname2 + '_gliderCTD', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()

    # ----------------------------------------------

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')


    #    # Plot the glider profiles
    scat3Dglider1 = ax.scatter(longlider1, latglider1, -depthglider1, s=5, c=chloroglider1, edgecolor='None',
                               norm=normchloro, cmap=cmapCTDchloro, zorder=4)

    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmap, orientation='vertical', pad=0.05, aspect=15, shrink=0.8, norm=normCTD,
                        extend='both')
    cbar.set_label('$mg\,m^{-3}$', fontname='Times New Roman', fontsize=18)
    cbar.set_clim(chloroCTDmin, chloroCTDmax)
    boundsCTDchloro = np.array((0.01, 0.03, 0.1, 0.3, 1.0))
    cbar.set_ticks(boundsCTDchloro)
    cbar.ax.set_yticklabels(boundsCTDchloro, fontname='Times New Roman', fontsize=16)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname2 + '_glider1', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()

    # ----------------------------------

    fig = plt.figure(num=None, figsize=(14, 8))
    ax = fig.gca(projection='3d')
    fig.patch.set_facecolor('white')

    # Plot the glider profiles
    scat3Dglider2 = ax.scatter(longlider2, latglider2, -depthglider2, s=5, c=chloroglider2, edgecolor='None',
                               norm=normchloro, cmap=cmapCTDchloro, zorder=4)

    # Add the colorbar
    cbar = fig.colorbar(scat3D, cmap=cmap, orientation='vertical', pad=0.05, aspect=15, shrink=0.8, norm=normCTD,
                        extend='both')
    cbar.set_label('$mg\,m^{-3}$', fontname='Times New Roman', fontsize=18)
    cbar.set_clim(chloroCTDmin, chloroCTDmax)
    boundsCTDchloro = np.array((0.01, 0.03, 0.1, 0.3, 1.0))
    cbar.set_ticks(boundsCTDchloro)
    cbar.ax.set_yticklabels(boundsCTDchloro, fontname='Times New Roman', fontsize=16)

    # Plot coastline
    # ax.plot3D(loncoast,latcoast,0,color='k',lw=0.5,zorder=6)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    # Change wall properties
    change_wall_prop(ax, coordinates, depths, angles)

    plt.savefig(figdir + figname2 + '_glider2', dpi=300, facecolor='w', edgecolor='w', transparent=False,
                bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()
