# Standard libs
import argparse
import os
import tempfile
from datetime import datetime
from glob import glob
from shutil import rmtree

# 3rd party libs
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pyart
import utm

from nexradaws.scripts.aws import get_nexrad_data
from utils import *

# Sweep to take
SWEEP = 0


def process_files(start_date, end_date, site, data_dir, out_dir):
    # Set up the directory for the radar download
    if data_dir is None:
        data_dir = tempfile.mkdtemp()
        tmp_dir = True
    else:
        tmp_dir = False

    # Start by trying to download the data (should go quick if the data is already in the right directory)
    data_dirs = get_nexrad_data(site, start_date, end_date, data_dir)

    # Make a list of the files to process
    files = []
    for dir in data_dirs:
        dir = os.path.join(dir, '*')
        files = files + glob(dir)

    # Do some bookkeeping for later
    image_base_dir = os.path.join(out_dir, 'images')
    nc_base_dir = os.path.join(out_dir, 'netcdf')

    # Iterate through the files
    first = True
    for f in files:
        print(f)
        radar = pyart.io.read_nexrad_archive(f)
        dt = pyart.graph.common.generate_radar_time_begin(radar)  # Scan time
        slice = radar.get_slice(SWEEP)

        if dt > datetime.strptime(end_date, "%Y%m%d-%H%M%S"):
            break

        if first:
            radar_lat = radar.latitude['data'][0]
            radar_lon = radar.longitude['data'][0]

            # Convert lat-lons to utm for
            x_radar, y_radar, _, _ = utm.from_latlon(radar_lat, radar_lon)

            # Init arrays for averaging
            phi_dp_running = np.zeros_like(radar.fields['differential_phase']['data'][slice].data)
            phi_dp_weighted_running = np.zeros_like(phi_dp_running)
            ref_linear_running = np.zeros_like(radar.fields['reflectivity']['data'][slice])
            eta_linear_running = np.zeros_like(ref_linear_running)

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
        ref.mask = gate_filter.gate_excluded[az_p, :]

        # Add data to the running sums
        phi_dp_running[~phi_dp.mask] += phi_dp[~phi_dp.mask]
        phi_dp_weighted_running[~phi_dp.mask] += phi_dp[~phi_dp.mask] * ref[~phi_dp.mask]
        ref_linear_running[~ref.mask] += db2pow(ref[~ref.mask])
        eta_linear_running[~ref.mask] += db2pow(ref[~ref.mask] + 11.6)

        # Make some plots
        plt.figure(figsize=(16, 8))
        plt.subplot(1, 2, 1)
        plt.xlim(-100, 100)
        plt.ylim(-100, 100)
        plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, phi_dp, vmin=0, vmax=360, cmap='nipy_spectral')
        plt.colorbar()
        plt.scatter(0, 0, color='k')

        plt.subplot(1, 2, 2)
        plt.xlim(-100, 100)
        plt.ylim(-100, 100)
        plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, ref, vmin=0, vmax=50)
        plt.colorbar()
        plt.scatter(0, 0, color='k')

        # Create directory if needed
        image_dir = os.path.join(image_base_dir, dt.strftime('%Y/%m/'))
        if not os.path.exists(image_dir): os.makedirs(image_dir)
        image_name = os.path.join(image_dir, dt.strftime('{site}_%Y%m%d_%H%M%S.png'.format(site=site)))

        plt.suptitle(dt.strftime('{site} %Y%m%d-%H%M%S Elev: {elev}'.format(elev=np.rad2deg(elev), site=site)))
        plt.savefig(image_name)
        plt.close()
        # plt.show(block=True)

    # Write out the netcdf
    start_time = datetime.strptime(start_date, "%Y%m%d-%H%M%S")
    if not os.path.exists(nc_base_dir): os.makedirs(nc_base_dir)
    nc_name = os.path.join(nc_base_dir, start_time.strftime("average_%Y%m%d.nc"))
    nc = netCDF4.Dataset(nc_name, 'w')

    # Add the dimensions
    az = nc.createDimension('az', size=(len(az_rad[:,0])),)
    rng = nc.createDimension('rng', size=(range_m.shape[1]))

    # Add the attributes
    attrs = {'num_scans': len(files),
             'start_time': start_date,
             'end_time': end_date,
             'radar_lat': radar_lat,
             'radar_lon': radar_lon,
             'sweep_number': SWEEP,
             'elevation': np.rad2deg(elev)
             }
    nc.setncatts(attrs)

    # Add the variables
    var = nc.createVariable('phi_dp_sum', datatype='f8', dimensions=('az', 'rng'))
    var.setncattr('units', 'degrees')
    var[:] = phi_dp_running

    var = nc.createVariable('phi_dp_weighted_sum', datatype='f8', dimensions=('az', 'rng'))
    var[:] = phi_dp_weighted_running

    var = nc.createVariable('ref_linear_sum', datatype='f8', dimensions=('az', 'rng'))
    var.setncattr('units', 'mm^6/m^3')
    var[:] = ref_linear_running

    var = nc.createVariable('eta_linear_sum', datatype='f8', dimensions=('az', 'rng'))
    var[:] = eta_linear_running

    var = nc.createVariable('azimuth', datatype='f8', dimensions=('az',))
    var.setncattr('units', 'radians')
    print(az_rad[az_p, 0])
    var[:] = az_rad[az_p, 0]

    var = nc.createVariable('range', datatype='f8', dimensions=('rng',))
    var.setncattr('units', 'm')
    var[:] = radar.range['data'][:]

    nc.close()

    # Delete the temporary folder if used
    if tmp_dir is True:
        rmtree(data_dir)

if __name__=='__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='start_date', help="YYYYmmdd-HHMMSS")
    parser.add_argument('-e', dest='end_date', help="YYYYmmdd-HHMMSS")
    parser.add_argument('-r', dest='radar')
    parser.add_argument('-d', dest='data_dir',
                        help='directory to download radar data (uses tmp dir otherwise)',
                        default=None)
    parser.add_argument('-o', dest='out_dir',
                        help='Directory to put images and netcdfs. Code will organize the dir structure')
    args = parser.parse_args()

    process_files(args.start_date, args.end_date, args.radar, args.data_dir, args.out_dir)



