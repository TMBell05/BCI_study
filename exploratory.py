import os
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pyart
import utm

from random_functions import save_data_as_image

# Sweep to process
sweep = 0

# Lat lon of Frio cave
lat_bat = 29.434507
lon_bat = -99.684534

files = glob('/Users/tbupper90/Data/nexrad/KDFX/raw/20160617/KDFX20160617_02*')

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
        phi_dp_running = np.zeros_like(radar.fields['differential_phase']['data'][slice].data)
        ref = np.zeros_like(radar.fields['reflectivity']['data'][slice])

        first = False

    # Get azimuth, range, and elevation for conversion to x and y
    range_m, az_deg = np.meshgrid(radar.range['data'], radar.azimuth['data'][slice])
    az_rad = np.deg2rad(az_deg)
    elev = np.deg2rad(np.mean(radar.elevation['data'][slice]))

    x_m = range_m * np.cos(elev) * np.sin(az_rad)
    y_m = range_m * np.cos(elev) * np.cos(az_rad)

    # Get the correct order of the azimuths
    az_p = np.argsort(az_rad, axis=0)[:,0]

    # Sort the x and y grids based on the order of the azimuths
    x_m = x_m[az_p, :]
    y_m = y_m[az_p, :]

    # # Apply filters and corrections to data
    gate_filter = pyart.filters.GateFilter(radar)
    gate_filter.exclude_below('differential_phase', 70.)
    gate_filter = pyart.correct.despeckle_field(radar, 'differential_phase', gatefilter=gate_filter)

    # Extract the desired data and get it in the correct order
    phi_dp = radar.fields['differential_phase']['data'][slice][az_p, :]
    ref = radar.fields['reflectivity']['data'][slice][az_p, :]

    # Convert the filter to mask
    phi_dp.mask = gate_filter.gate_excluded[az_p, :]

    # Get differential phase data
    phi_dp_running[~phi_dp.mask] += phi_dp[~phi_dp.mask] * ref[~phi_dp.mask]

    # plt.figure(figsize=(16, 8))
    # plt.subplot(1, 2, 1)
    # plt.xlim(-80, 80)
    # plt.ylim(-80, 80)
    # plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, phi_dp, vmin=0, vmax=360, cmap='nipy_spectral')
    # plt.colorbar()
    # plt.scatter(0, 0, color='k')
    #
    # plt.subplot(1, 2, 2)
    # plt.xlim(-80, 80)
    # plt.ylim(-80, 80)
    # plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, ref, vmin=0, vmax=50)
    # plt.colorbar()
    # plt.scatter(0, 0, color='k')
    #
    # plt.suptitle(start_time.strftime('KDFX %Y%m%d-%H%M%S Elev: {elev}'.format(elev=np.rad2deg(elev))))
    # plt.savefig(f.split('/')[-1] + '.png')
    # plt.close()
    # # plt.show(block=True)

count = len(files)

# x_m -= x_rel_m
# y_m -= y_rel_m

# create quick plot
plt.figure()
plt.xlim(-100, 100)
plt.ylim(-100, 100)
plt.pcolormesh(x_m*1e-3, y_m*1e-3, phi_dp_running/count, vmin=0, vmax=360, cmap='nipy_spectral')
plt.colorbar()
plt.scatter(0, 0, color='y')
# plt.title(start_time.strftime('KDFX %Y%m%d-%H%M%S Elev: {elev}'.format(elev=np.rad2deg(elev))))
plt.savefig(start_time.strftime('avg_KDFX_%Y%m%d.png'))
# plt.show()
plt.close()


save_data_as_image(x_m*1e-3, y_m*1e-3, phi_dp_running/count, 'test2.png')
