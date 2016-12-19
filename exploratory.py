import os
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pyart
import utm

# Sweep to process
sweep = 0

# Lat lon of Frio cave
lat_bat = 29.434507
lon_bat = -99.684534

files = glob('/Users/tylerbell/data/nexrad/KDFX/raw/20160617/KDFX20160617_02*')

first = True
for f in files:
    print f
    radar = pyart.io.read_nexrad_archive(f)
    start_time = pyart.graph.common.generate_radar_time_begin(radar)
    slice = radar.get_slice(sweep)

    if first:
        radar_lat = radar.latitude['data'][0]
        radar_lon = radar.longitude['data'][0]

        # Convert lat-lons to utm for
        x_bat, y_bat, _, _ = utm.from_latlon(lat_bat, lon_bat)
        x_radar, y_radar, _, _ = utm.from_latlon(radar_lat, radar_lon)

        # Calc relative x and y of roost
        x_rel_m = (x_bat - x_radar)
        y_rel_m = (y_bat - y_radar)

        # Get rid of misc crap
        del x_bat, y_bat, x_radar, y_radar, _

        # Init arrays for averaging
        kdp_hv = np.zeros_like(radar.fields['differential_phase']['data'][slice].data)
        ref = np.zeros_like(radar.fields['reflectivity']['data'][slice])

        first = False

    # Get azimuth, range, and elevation for conversion to x and y
    range_m, az_deg = np.meshgrid(radar.range['data'], radar.azimuth['data'][slice])
    az_rad = np.deg2rad(az_deg)
    elev = np.deg2rad(np.mean(radar.elevation['data'][slice]))

    x_m = range_m * np.cos(elev) * np.sin(az_rad)
    y_m = range_m * np.cos(elev) * np.cos(az_rad)

    tmp = radar.fields['differential_phase']['data'][slice]
    ind = tmp.data > 70


    # Get differential phase data
    kdp_hv[ind] += tmp[ind]
    ref += radar.fields['reflectivity']['data'][slice]

    plt.figure()
    plt.xlim(-80, 80)
    plt.ylim(-80, 80)
    plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, kdp_hv, vmin=0, vmax=360, cmap='nipy_spectral')
    plt.colorbar()
    plt.scatter(0, 0, color='k')
    # plt.title(start_time.strftime('KDFX %Y%m%d-%H%M%S Elev: {elev}'.format(elev=np.rad2deg(elev))))
    # plt.savefig(test_file.split('/')[-1] + '.png')
    plt.show()

count = len(files)

# x_m -= x_rel_m
# y_m -= y_rel_m

# create quick plot
plt.figure()
plt.xlim(-40, 40)
plt.ylim(-40, 40)
plt.pcolormesh(x_m*1e-3, y_m*1e-3, kdp_hv, vmin=0, vmax=360, cmap='nipy_spectral')
plt.colorbar()
plt.scatter(0, 0, color='y')
# plt.title(start_time.strftime('KDFX %Y%m%d-%H%M%S Elev: {elev}'.format(elev=np.rad2deg(elev))))
# plt.savefig(test_file.split('/')[-1] + '.png')
plt.show()