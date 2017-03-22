# Standard libs
import argparse
import logging
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

# Cave file
cave_csv = 'cave_locations.csv'
caves = read_caves(cave_csv, 'KDFX')

# Expected radar grid (azimuth and range)
AZ_SIZE = None
RNG_SIZE = None
AZ = None
RNG = None


def process_files(start_date, end_date, site, data_dir, out_dir, verbose=False):
    if verbose:
        logging.basicConfig(format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG)
    else:
        # Set up logging
        logging.basicConfig(format="%(asctime)s:%(levelname)s:%(message)s", level=logging.INFO)

    # Set up the directory for the radar download
    if data_dir is None:
        data_dir = tempfile.mkdtemp()
        logging.debug("Created temp dir: {0}".format(data_dir))
        tmp_dir = True
    else:
        tmp_dir = False

    # Start by trying to download the data (should go quick if the data is already in the right directory)
    logging.debug("Starting to download data")
    data_dirs = get_nexrad_data(site, start_date, end_date, data_dir)
    # Make a list of the files to process
    files = []
    for dir in data_dirs:
        dir = os.path.join(dir, '*')
        files = files + glob(dir)
    logging.debug("Processing {0} files".format(len(files)))

    # Do some bookkeeping for later
    image_base_dir = os.path.join(out_dir, 'images')
    logging.debug("Images being written to {0}".format(image_base_dir))
    nc_base_dir = os.path.join(out_dir, 'netcdf')
    logging.debug("NetCDFs being written to {0}".format(nc_base_dir))

    # Iterate through the files
    first = True
    count = 0
    for f in files:
        logging.info("Processing {0}".format(f))

        try:
            radar = pyart.io.read_nexrad_archive(f)
        except Exception as e:
            print("Can't open file " + f)

        dt = pyart.graph.common.generate_radar_time_begin(radar)  # Scan time
        slice = radar.get_slice(SWEEP)

        if radar.metadata['vcp_pattern'] not in [31, 32]:
            logging.warning("VCP other than clear air mode found, skipping file")
            continue

        if dt > datetime.strptime(end_date, "%Y%m%d-%H%M%S"):
            break

        if first:
            radar_lat = radar.latitude['data'][0]
            radar_lon = radar.longitude['data'][0]

            cave_x = []
            cave_y = []
            # Convert cave lat/lon to x/y coord system
            x_radar, y_radar, _, _ = utm.from_latlon(radar_lat, radar_lon)

            for cave in caves:
                # Convert lat-lons to utm for
                x_bat, y_bat, _, _ = utm.from_latlon(float(cave['lat']), float(cave['lon']))
                # Calc relative x and y of roost
                x_rel_m = (x_bat - x_radar)
                y_rel_m = (y_bat - y_radar)
                cave['x'] = x_rel_m / 1e3
                cave['y'] = y_rel_m / 1e3

                cave_x.append(x_rel_m/1e3)
                cave_y.append(y_rel_m/1e3)
                # Get rid of misc crap
                del x_bat, y_bat, _

            # Convert lat-lons to utm for
            x_radar, y_radar, _, _ = utm.from_latlon(radar_lat, radar_lon)

            # Init arrays for averaging
            phi_dp_running = np.zeros_like(radar.fields['differential_phase']['data'][slice].data)
            phi_dp_weighted_running = np.zeros_like(phi_dp_running)
            ref_linear_running = np.zeros_like(radar.fields['reflectivity']['data'][slice])
            eta_linear_running = np.zeros_like(ref_linear_running)

        # Get azimuth, range, and elevation for conversion to x and y
        range_m, az_deg = np.meshgrid(radar.range['data'], radar.azimuth['data'][slice])
        az_rad = np.deg2rad(az_deg)
        elev = np.deg2rad(np.mean(radar.elevation['data'][slice]))

        x_m = range_m * np.cos(elev) * np.sin(az_rad)
        y_m = range_m * np.cos(elev) * np.cos(az_rad)

        # Get the correct order of the azimuths
        az_p = np.argsort(az_rad, axis=0)[:, 0]

        # Sort the x and y grids based on the order of the azimuths
        x_m = x_m[az_p, :]
        y_m = y_m[az_p, :]

        if first:
            RNG = range_m[0, :]
            RNG_SIZE = RNG.size
            AZ = az_rad[:, 0][az_p]
            AZ_SIZE = AZ.size
            first = False

        # # Apply filters and corrections to data
        logging.debug("Applying corrections")
        gate_filter = pyart.filters.GateFilter(radar)
        gate_filter.exclude_below('differential_phase', 70.)
        gate_filter = pyart.correct.despeckle_field(radar, 'differential_phase', gatefilter=gate_filter)

        # Extract the desired data and get it in the correct order
        phi_dp = radar.fields['differential_phase']['data'][slice][az_p, :]
        ref = radar.fields['reflectivity']['data'][slice][az_p, :]

        # Convert the filter to mask
        phi_dp.mask = gate_filter.gate_excluded[az_p, :]
        ref.mask = gate_filter.gate_excluded[az_p, :]

        # Check to make sure the data lines up with the desired size
        if phi_dp_running.shape != phi_dp.shape:
            logging.warning("Data size other than ({}, {}) found".format(AZ_SIZE, RNG_SIZE))
            logging.warning("Data size: ({}, {})".format(phi_dp.shape[0], phi_dp.shape[1]))
            if az_rad.size < AZ_SIZE:
                logging.debug("Applying pad to account for having too few azimuths")
                diff = abs(az_rad.size - AZ_SIZE)
                phi_dp = np.pad(phi_dp, diff, mode='constant')[diff:, diff:-diff]
                ref = np.pad(ref, diff, mode='constant')[diff:, diff:-diff]
            elif az_rad.size > AZ_SIZE:
                logging.debug("Chopping off end of grid because too many azimuths")
                phi_dp = phi_dp[:AZ_SIZE, :]
                ref = ref[:AZ_SIZE, :]
            else:
                logging.critical("SOMETHING WENT WRONG WITH THE RANGE")
                raise Exception

        # Add data to the running sums
        # try:
        logging.debug("Added data to running sums")
        phi_dp_running[~phi_dp.mask] += phi_dp[~phi_dp.mask]
        phi_dp_weighted_running[~phi_dp.mask] += phi_dp[~phi_dp.mask] * ref[~phi_dp.mask]
        ref_linear_running[~ref.mask] += db2pow(ref[~ref.mask])
        eta_linear_running[~ref.mask] += db2pow(ref[~ref.mask] + 11.6)
        # except IndexError as e:
        #     # TODO- account for this somehow other than an exception
        #     print("Index error Encountered, taking only first 360 lines on axis zero")
        #     phi_dp_running[~phi_dp.mask[0:359, :]] += phi_dp[~phi_dp.mask[0:359, :]]
        #     phi_dp_weighted_running[~phi_dp.mask[0:359, :]] += phi_dp[~phi_dp.mask[0:359, :]] * ref[~phi_dp.mask[0:359, :]]
        #     ref_linear_running[~ref.mask[0:359, :]] += db2pow(ref[~ref.mask[0:359, :]])
        #     eta_linear_running[~ref.mask[0:359, :]] += db2pow(ref[~ref.mask[0:359, :]] + 11.6)
        #
        #     az_p = az_p[0:359]
        # except Exception as e:
        #     print("Unknown error: " + e)
        #     continue

        # Make some plots
        plt.figure(figsize=(16, 8))
        plt.subplot(1, 2, 1)
        plt.xlim(-100, 100)
        plt.ylim(-100, 100)
        plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, phi_dp, vmin=0, vmax=360, cmap='nipy_spectral')
        plt.colorbar()
        plt.scatter(0, 0, color='k')
        plt.scatter(cave_x, cave_y, c='y')

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
        logging.debug("Saving image: {}".format(image_name))
        plt.savefig(image_name)
        plt.close()
        # plt.show(block=True)

        count += 1

    # If no files were processed
    if count == 0:
        return

    # Write out the netcdf
    logging.info("Preparing netCDF")
    start_time = datetime.strptime(start_date, "%Y%m%d-%H%M%S")
    if not os.path.exists(nc_base_dir): os.makedirs(nc_base_dir)
    nc_name = os.path.join(nc_base_dir, start_time.strftime("{}_average_%Y%m%d.nc".format(site)))
    nc = netCDF4.Dataset(nc_name, 'w')

    # Add the dimensions
    print(phi_dp_running.shape[1], AZ.shape)
    az = nc.createDimension('az', size=phi_dp_running.shape[0],)
    rng = nc.createDimension('rng', size=phi_dp_running.shape[1],)


    # Add the attributes
    attrs = {'num_scans': count,
             'start_time': start_date,
             'end_time': end_date,
             'radar_lat': radar_lat,
             'radar_lon': radar_lon,
             'sweep_number': SWEEP,
             'elevation': np.rad2deg(elev),
             'site': site
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
    var[:] = AZ

    var = nc.createVariable('range', datatype='f8', dimensions=('rng',))
    var.setncattr('units', 'm')
    var[:] = RNG

    nc.close()
    logging.info("NetCDF write successful")

    # Delete the temporary folder if used
    if tmp_dir is True:
        logging.debug("Removing temp dir")
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



