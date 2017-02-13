import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.interpolate import splprep, splev

import medgib_functions

# Files and directories

datadir = '/home/ctroupin/Publis/20150528_MedGib/data/'
datafiles = sorted(glob.glob(datadir + '*csv'))
coastdir = '/home/ctroupin/DataOceano/Coastlines/'
coastfile = 'medsea_coast2b.dat'

currentdir = '/home/ctroupin/Publis/20150528_MedGib/python/'
currentfile1 = 'AtlanticJet2.dat'

figdir = '/home/ctroupin/public_html/MedGib/'
valex = -999.

plotlarge, plotmedium, plotsmall = 0, 1, 0

coordinates1 = np.array((-7, 9., 34.5, 40))
coordinates2 = np.array((-7, 2., 34.5, 37.5))
coordinates3 = np.array((-7, -3., 35, 36.75))
dlon1, dlat1 = 2.0, 1.0
dlon2, dlat2 = 2.0, 1.0
dlon3, dlat3 = 1.0, .5

res1, res2, res3 = 'l', 'l', 'l'

sstmin, sstmax = 16., 25.
cmapsst = plt.cm.RdYlBu_r
normsst = colors.Normalize(vmin=sstmin, vmax=sstmax)
boundsst = np.arange(sstmin, sstmax + .1, 1.0)


# Load and smoorth currents
lonJet, latJet = np.loadtxt(currentdir + currentfile1, usecols=(0, 1), unpack=True)
t = np.linspace(0, 1, len(lonJet))
t2 = np.linspace(0, 1, 4 * len(lonJet))

# spline parameters
s = 0.01  # smoothness parameter
k = 3  # spline order
nest = 4  # estimate of number of knots needed (-1 = maximal)

# find the knot points
tckp, u = splprep([t, lonJet, latJet], s=s, k=k, nest=-1)

# evaluate spline, including interpolated points
xnew, lonJetsmooth, latJetsmooth = splev(t2, tckp)

if plotlarge:

    # Figure 1 - Large domain
    # -----------------------

    figname = 'Medgib_SST_large_domain'

    # Prepare map and projection
    fig, m, ax = medgib_functions.prepare_map(coordinates1, res1)

    # Load coast
    lonc, latc = medgib_functions.medgib_load_coast(coastdir, coastfile, valex)
    lonc, latc = m(lonc, latc)
    m.plot(lonc, latc, 'k-', lw=0.5, zorder=4)

    # Plot currents
    lonJetsmooth, latJetsmooth = m(lonJetsmooth, latJetsmooth)
    plt.plot(lonJetsmooth, latJetsmooth, 'k-', lw=2, zorder=5)

    # Add SST measured by buoys
    for f in datafiles:
        # print f
        file2load = f
        code, year, month, day, hour, minute, second, lat, lon, sst,\
        u, v, lat_qc, lon_qc, sst_qc, V_QC = np.loadtxt(file2load, delimiter=',', unpack=True)

        # Keep only good positions
        lon[lon_qc > 1] = np.nan
        lon, lat = m(lon, lat)
        sstscat = plt.scatter(lon, lat, s=5, c=sst, edgecolor='none', cmap=cmapsst, norm=normsst)

    if len(datafiles) > 0:
        # colorbar
        cbar = fig.colorbar(sstscat, cmap=cmapsst, norm=normsst,
                            orientation='horizontal', pad=0.025,
                            aspect=20, shrink=1, extend='both')
        cbar.set_ticks(boundsst)
        cbar.set_label(r'^{\circ}C', fontsize=16, rotation=0)

    m.drawparallels(np.arange(np.round_(coordinates1[2]), coordinates1[3], dlat1), linewidth=0.5,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(coordinates1[0], coordinates1[1], dlon1), linewidth=0.5,
                    labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmapscale(-4., 38.5, -4., 38, 200, barstyle='simple', units='km', fontsize=12, zorder=4)

    m.fillcontinents(ax=ax, color='w', zorder=3)
    # plt.savefig(figdir+figname, dpi=300, facecolor='w', edgecolor='w',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close()

if plotmedium:

    # Figure 2 - Medium domain (14 days)
    # -----------------------

    figname = 'Medgib_SST_medium_domain_currents3'

    # Prepare map and projection
    fig, m, ax = medgib_functions.prepare_map(coordinates2, res2)

    # Load coast
    lonc, latc = medgib_functions.medgib_load_coast(coastdir, coastfile, valex)
    lonc, latc = m(lonc, latc)
    m.plot(lonc, latc, 'k-', lw=0.5, zorder=4)

    # Plot currents
    lonJetsmooth, latJetsmooth = m(lonJetsmooth, latJetsmooth)
    plt.plot(lonJetsmooth, latJetsmooth, 'k-', lw=2, zorder=5)

    # Define radii for the gyres
    WAGradius, EAGradius, Smallradius = 55., 49., 18.
    WAGlon, WAGlat = -4.25, 35.75
    EAGlon, EAGlat = -2.1, 35.6
    Smalllon, Smalllat = -3.25, 35.6
    medgib_functions.equi(m, WAGlon, WAGlat, WAGradius, lw=2., color='k', alpha=0.7)  # WAG
    medgib_functions.equi(m, EAGlon, EAGlat, EAGradius, lw=2., color='k', alpha=0.7)  # EAG
    medgib_functions.equi(m, Smalllon, Smalllat, Smallradius, lw=2., color='k', alpha=0.7)

    # Add text
    plt.text(245500, 138862, 'WAG', va='center', ha='center')
    plt.text(453303, 129610, 'EAG', va='center', ha='center')
    plt.text(302973, 215000, 'AJ', va='center', ha='center')

    # Add arrows
    # ----------

    myarrowdict = dict(arrowstyle="-|>", fc="k")

    # WAG
    glon1, glat1, baz = medgib_functions.shoot(WAGlon, WAGlat, -89, WAGradius)
    glon2, glat2, baz = medgib_functions.shoot(WAGlon, WAGlat, -90, WAGradius)
    x, y = m(glon1, glat1)
    x2, y2 = m(glon2, glat2)
    plt.annotate('', xy=(x, y), xytext=(x2, y2), arrowprops=myarrowdict)

    glon1, glat1, baz = medgib_functions.shoot(WAGlon, WAGlat, 92, WAGradius)
    glon2, glat2, baz = medgib_functions.shoot(WAGlon, WAGlat, 90, WAGradius)
    x, y = m(glon1, glat1)
    x2, y2 = m(glon2, glat2)
    plt.annotate('', xy=(x, y), xytext=(x2, y2), arrowprops=myarrowdict)

    # EAG
    glon1, glat1, baz = medgib_functions.shoot(EAGlon, EAGlat, -89, EAGradius)
    glon2, glat2, baz = medgib_functions.shoot(EAGlon, EAGlat, -90, EAGradius)
    x, y = m(glon1, glat1)
    x2, y2 = m(glon2, glat2)
    plt.annotate('', xy=(x, y), xytext=(x2, y2), arrowprops=myarrowdict)

    glon1, glat1, baz = medgib_functions.shoot(EAGlon, EAGlat, 92, EAGradius)
    glon2, glat2, baz = medgib_functions.shoot(EAGlon, EAGlat, 90, EAGradius)
    x, y = m(glon1, glat1)
    x2, y2 = m(glon2, glat2)
    plt.annotate('', xy=(x, y), xytext=(x2, y2), arrowprops=myarrowdict)

    # Small eddy
    glon1, glat1, baz = medgib_functions.shoot(Smalllon, Smalllat, -90, Smallradius)
    glon2, glat2, baz = medgib_functions.shoot(Smalllon, Smalllat, -89, Smallradius)
    x, y = m(glon1, glat1)
    x2, y2 = m(glon2, glat2)
    plt.annotate('', xy=(x, y), xytext=(x2, y2), arrowprops=myarrowdict)

    glon1, glat1, baz = medgib_functions.shoot(Smalllon, Smalllat, 89, Smallradius)
    glon2, glat2, baz = medgib_functions.shoot(Smalllon, Smalllat, 90, Smallradius)
    x, y = m(glon1, glat1)
    x2, y2 = m(glon2, glat2)
    plt.annotate('', xy=(x, y + 1000.), xytext=(x2, y2 + 1000.), arrowprops=myarrowdict)

    # Atlantic Jet
    plt.annotate('', xy=(255000, 203052), xytext=(250000, 203052), arrowprops=myarrowdict)
    plt.annotate('', xy=(340000, 90000), xytext=(335000, 90000), arrowprops=myarrowdict)
    plt.annotate('', xy=(444000, 178000), xytext=(440000, 178000), arrowprops=myarrowdict)
    plt.annotate('', xy=(532000, 137500), xytext=(530000, 137500), arrowprops=myarrowdict)

    # Add SST measured by buoys
    for f in datafiles:
        # print f
        file2load = f
        code, year, month, day, hour, minute, second, lat, lon, sst, u, v, lat_qc, lon_qc, sst_qc, V_QC = np.loadtxt(
            file2load, delimiter=',', unpack=True)

        goodday = np.where((day <= 25) & (month == 9))[0]
        lon, lat, sst, lon_qc = lon[goodday], lat[goodday], sst[goodday], lon_qc[goodday]

        # Keep only good positions
        lon[lon_qc > 1] = np.nan
        lon, lat = m(lon, lat)
        sstscat = plt.scatter(lon, lat, s=5, c=sst, edgecolor='none', cmap=cmapsst, norm=normsst)

    # colorbar
    cbar = fig.colorbar(sstscat, cmap=cmapsst, norm=normsst, orientation='horizontal', pad=0.025, aspect=20, shrink=1,
                        extend='both')
    cbar.set_ticks(boundsst)
    cbar.set_label(r'^{\circ}C', fontsize=16, rotation=0)

    m.drawparallels(np.arange(np.round_(coordinates2[2]), coordinates2[3], dlat2), linewidth=0.5,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(coordinates2[0], coordinates2[1], dlon2), linewidth=0.5,
                    labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)

    m.fillcontinents(ax=ax, color='w', zorder=3)
    plt.savefig(figdir + figname, dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close()

if plotsmall:

    # Figure 3 - Small domain (4 days)
    # -----------------------

    figname = 'Medgib_SST_small_domain'

    # Prepare map and projection
    fig, m, ax = medgib_functions.prepare_map(coordinates3, res3)

    # Load coast
    lonc, latc = medgib_functions.medgib_load_coast(coastdir, coastfile, valex)
    lonc, latc = m(lonc, latc)
    m.plot(lonc, latc, 'k-', lw=0.5, zorder=4)

    # Add SST measured by buoys
    for f in datafiles:
        # print f
        file2load = f
        code, year, month, day, hour, minute, second, lat, lon, sst, u, v, lat_qc, lon_qc, sst_qc, V_QC = np.loadtxt(
            file2load, delimiter=',', unpack=True)

        goodday = np.where((day <= 13) & (month == 9))[0]
        lon, lat, sst, lon_qc = lon[goodday], lat[goodday], sst[goodday], lon_qc[goodday]

        # Keep only good positions
        lon[lon_qc > 1] = np.nan
        lon, lat = m(lon, lat)
        sstscat = plt.scatter(lon, lat, s=5, c=sst, edgecolor='none', cmap=cmapsst, norm=normsst)

    # colorbar
    cbar = fig.colorbar(sstscat, cmap=cmapsst, norm=normsst, orientation='horizontal', pad=0.025, aspect=20, shrink=1,
                        extend='both')
    cbar.set_ticks(boundsst)
    cbar.set_label(r'^{\circ}C', fontsize=16, rotation=0)

    m.drawparallels(np.arange(np.round_(coordinates2[2]), coordinates2[3], dlat3), linewidth=0.5,
                    labels=[1, 0, 0, 0], fontname='Times New Roman', fontsize=16, zorder=1)
    m.drawmeridians(np.arange(coordinates2[0], coordinates2[1], dlon3), linewidth=0.5,
                    labels=[0, 0, 1, 0], fontname='Times New Roman', fontsize=16, zorder=1)

    m.fillcontinents(ax=ax, color='w', zorder=3)
    plt.savefig(os.path.join(figdir, figname), dpi=300, facecolor='w', edgecolor='w',
                transparent=False, bbox_inches='tight', pad_inches=0.1)
    # plt.show()
    plt.close()
