import os
import logging
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from scipy import interpolate

from utils import pow2db

# Set up logging
logging.basicConfig(format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG)

# Specify the seasons desired
seasons = [2013, 2014, 2015, 2016]
radars = ['KDFX', 'KEWX', 'KGRK', 'KSJT']

# Specify the RCS of the bat
RCS = 11.6  # cm^2

# Levels
levels = [x for x in range(0, 301, 50)]
print(levels)
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
            tmp_ref = np.ma.masked_where(nc['ref_linear_sum'][:] < contour_value, nc['ref_linear_sum'])
            tmp_eta = np.ma.MaskedArray(nc['eta_linear_sum'], mask=tmp_ref.mask)

            keepers = nc['ref_linear_sum'][:] > contour_value

            # plt.figure()
            # plt.pcolormesh(x_m * 1.e-3, y_m * 1.e-3, pow2db(tmp_ref), vmin=0, vmax=50)
            # plt.colorbar()
            # plt.xlim(-200, 200)
            # plt.ylim(-200, 200)
            # plt.savefig('data/season_avgs/raw/{}_{}.png'.format(radar, nc.start_time))
            # plt.close()

            linear_ref += tmp_ref.filled(0)
            linear_eta += tmp_eta.filled(0)

            nc.close()

        # Do some QC?
        linear_ref = np.ma.masked_where(linear_ref == 0, linear_ref)
        avg_ref = pow2db(linear_ref/float(len(files)))
        print(linear_eta)

        plt.figure()
        plt.pcolormesh(x_m / 1000., y_m / 1000., avg_ref, vmin=0, vmax=50)
        plt.xlim(-200, 200)
        plt.ylim(-200, 200)
        plt.colorbar()
        plt.title('{} {} Reflectivity'.format(radar, season))
        plt.savefig('data/season_avgs/ref_{}_{}_avg.png'.format(radar, season))
        plt.close()
        # plt.show()

        linear_eta = np.ma.masked_where(linear_eta == 0, linear_eta)
        plt.figure()
        plt.pcolormesh(x_m / 1000., y_m / 1000., linear_eta/float(len(files))/RCS, vmax=100)
        plt.xlim(-200, 200)
        plt.ylim(-200, 200)
        plt.colorbar()
        # plt.contour(x_m / 1000., y_m / 1000., linear_eta/float(len(files))/RCS, levels=levels, colors='k')
        plt.title('{} {} Nbats/km^3'.format(radar, season))
        plt.savefig('data/season_avgs/sigma_{}_{}_avg.png'.format(radar, season))
        plt.close()
        # plt.show()




