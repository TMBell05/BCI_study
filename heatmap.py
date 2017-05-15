import os
import logging
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

from utils import pow2db

# Set up logging
logging.basicConfig(format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG)

# Specify the seasons desired
seasons = [2013, 2014, 2015, 2016]
radars = ['KDFX', 'KEWX', 'KGRK', 'KSJT']

# Location of daily averages
source_dir = 'data/netcdf_daily_average'

# For testing purposes
# seasons = [seasons[-1]]
# radars = [radars[0]]

for season in seasons:
    for radar in radars:

        # Get the files
        search_str = os.path.join(source_dir, "{}*{}*.nc".format(radar, season))
        files = glob(search_str)

        # If here, this is the first round
        first = True
        for f in files:
            print(f)
            nc = netCDF4.Dataset(f)

            # Get the grid figured out
            range_m, az_rad = np.meshgrid(nc['range'][:], nc['azimuth'][:])
            elev = np.deg2rad(nc.elevation)

            x_m = range_m * np.cos(elev) * np.sin(az_rad)
            y_m = range_m * np.cos(elev) * np.cos(az_rad)

            if first:
                linear_ref = np.zeros_like(nc['ref_linear_sum'])
                linear_eta = np.zeros_like(nc['eta_linear_sum'])

            contour_value = np.max(pow2db(nc['ref_linear_sum'][:] / nc.num_scans))/7.5
            keepers = nc['ref_linear_sum'][:] > contour_value

            linear_ref[keepers] += nc['ref_linear_sum'][:][keepers]
            linear_eta[keepers] += nc['eta_linear_sum'][:][keepers]

            nc.close()

        linear_ref = np.ma.masked_where(linear_ref == 0, linear_ref)
        print(linear_ref)
        plt.figure()
        plt.pcolormesh(x_m / 1000., y_m / 1000., pow2db(linear_ref/float(len(files))))
        plt.xlim(-200, 200)
        plt.ylim(-200, 200)
        plt.colorbar()
        plt.savefig('data/season_avgs/{}_{}_avg.png'.format(radar, season))
        plt.close()
        # plt.show()




