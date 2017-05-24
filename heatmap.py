import os
import logging
from glob import glob

import cartopy.crs as ccrs
import gdal
from osr import SpatialReference
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import pyproj
import utm
from scipy import interpolate

from utils import pow2db, interpolate

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
        search_str = os.path.join(source_dir, "{}*{}*.nc".format(radar, season))
        files = glob(search_str)

        # If here, this is the first round
        first = True
        for f in files:
            try:
                print(f)
                nc = netCDF4.Dataset(f)

                # Get the grid figured out
                range_m, az_rad = np.meshgrid(nc['range'][:], nc['azimuth'][:])
                elev = np.deg2rad(nc.elevation)

                x_m = range_m * np.cos(elev) * np.sin(az_rad)
                y_m = range_m * np.cos(elev) * np.cos(az_rad)

                if first:
                    # Set up grid to interpolate to
                    XI = np.arange(-460000., 460000., 1000.)
                    YI = np.arange(460000., -460000., -1000.)
                    XI, YI = np.meshgrid(XI, YI)

                    # Create running sum arrays
                    linear_ref = np.zeros_like(XI)
                    linear_eta = np.zeros_like(XI)

                    # Get the AZ and RNG of the grid
                    AZI = np.arctan2(YI, XI)
                    AZI = -AZI + np.pi/2.
                    AZI = np.where(AZI < 0, AZI + 2.*np.pi, AZI)
                    RNGI = np.sqrt(XI**2. + YI**2.)

                    # Coord transform stufffff
                    x_radar, y_radar, zone, _ = utm.from_latlon(nc.radar_lat, nc.radar_lon)
                    utm_xi = x_radar + XI
                    utm_yi = y_radar + YI

                    # proj4 = '+proj=utm +zone={} +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs'.format(zone)
                    proj4 = '+ellps=WGS84 +proj=ortho +lon_0={} +lat_0={} +no_defs'.format(nc.radar_lon, nc.radar_lat)
                    proj = pyproj.Proj(proj4)

                    lon, lat = proj(XI, YI, inverse=True)

                    # Set to false so this doesn't repeat for each file
                    first = False

                contour_value = np.max(pow2db(nc['ref_linear_sum'][:] / nc.num_scans))/7.5
                tmp_ref = np.ma.masked_where(nc['ref_linear_sum'][:] < contour_value, nc['ref_linear_sum'])
                tmp_eta = np.ma.MaskedArray(nc['eta_linear_sum'], mask=tmp_ref.mask)

                sort = np.argsort(nc['azimuth'][:])

                new_ref = interpolate(nc['azimuth'][:], nc['range'][:], tmp_ref, AZI, RNGI)
                new_eta = interpolate(nc['azimuth'][:], nc['range'][:], tmp_eta, AZI, RNGI)

                keepers = nc['ref_linear_sum'][:] > contour_value

                plt.figure()
                # plt.pcolormesh(lon, lat, pow2db(new_eta), vmin=0, vmax=100)
                plt.pcolormesh(XI * 1.e-3, YI * 1.e-3, pow2db(new_eta), vmin=0, vmax=100)
                plt.colorbar()
                plt.xlim(-200, 200)
                plt.ylim(-200, 200)
                plt.savefig('data/season_avgs/raw/{}_{}.png'.format(radar, nc.start_time))
                # plt.show()
                plt.close()

                linear_ref += new_ref
                linear_eta += new_eta

                nc.close()
            except Exception as e:
                print(e)
                print("Exception for file {}".format(f))

        # Do some QC?
        linear_ref = np.ma.masked_where(linear_ref == 0, linear_ref)
        avg_ref = pow2db(linear_ref/float(len(files)))

        plt.figure()
        plt.pcolormesh(XI / 1000., YI / 1000., avg_ref, vmin=0, vmax=50)
        plt.xlim(-200, 200)
        plt.ylim(-200, 200)
        plt.colorbar()
        plt.title('{} {} Reflectivity'.format(radar, season))
        plt.savefig('data/season_avgs/ref_{}_{}_avg.png'.format(radar, season))
        plt.close()
        # plt.show()

        linear_eta = np.ma.masked_where(linear_eta == 0, linear_eta)
        avg_eta = linear_eta/float(len(files))/RCS
        plt.figure()
        plt.pcolormesh(XI / 1000., YI / 1000., avg_eta, vmax=500)
        plt.xlim(-200, 200)
        plt.ylim(-200, 200)
        plt.colorbar()
        # plt.contour(x_m / 1000., y_m / 1000., linear_eta/float(len(files))/RCS, levels=levels, colors='k')
        plt.title('{} {} Nbats/km^3'.format(radar, season))
        plt.savefig('data/season_avgs/sigma_{}_{}_avg.png'.format(radar, season))
        plt.close()
        # plt.show()

        # Output geotiff
        xmin, ymin, xmax, ymax = [np.min(XI), np.min(YI), np.max(XI), np.max(YI)]
        nx = avg_ref.shape[0]
        ny = avg_ref.shape[1]
        xres = 1000.
        yres = 1000.

        geotransform = (xmin, xres, 0, ymax, 0, -yres)

        dst_ds = gdal.GetDriverByName('GTiff').Create('data/season_avgs/sigma_{}_{}_avg.tiff'.format(radar, season),
                                                      ny, nx, 1, gdal.GDT_Float64)
        dst_ds.SetGeoTransform(geotransform)
        srs = SpatialReference()
        srs.ImportFromProj4(proj4)
        dst_ds.SetProjection(srs.ExportToWkt())
        dst_ds.GetRasterBand(1).WriteArray(avg_eta.filled(-9999))
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.FlushCache()

        dst_ds = gdal.GetDriverByName('GTiff').Create('data/season_avgs/ref_{}_{}_avg.tiff'.format(radar, season),
                                                      ny, nx, 1, gdal.GDT_Float64)
        dst_ds.SetGeoTransform(geotransform)
        srs = SpatialReference()
        srs.ImportFromProj4(proj4)
        dst_ds.SetProjection(srs.ExportToWkt())
        dst_ds.GetRasterBand(1).WriteArray(avg_ref.filled(-9999))
        dst_ds.GetRasterBand(1).SetNoDataValue(-9999)
        dst_ds.FlushCache()


