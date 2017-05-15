import csv
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import os
import utm
from shapely.geometry import Polygon
from glob import glob
from utils import greens_theorem, read_caves, pow2db, distance

MIN_AREA = 75  # km^2
MAX_AREA = 12000  # km^2
MAX_DISTANCE = 75  # km^2


def get_sba(nc_file, cave_csv):
    nc = netCDF4.Dataset(nc_file)

    # Get the grid figured out
    range_m, az_rad = np.meshgrid(nc['range'][:], nc['azimuth'][:])
    elev = np.deg2rad(nc.elevation)

    x_m = range_m * np.cos(elev) * np.sin(az_rad)
    y_m = range_m * np.cos(elev) * np.cos(az_rad)

    plt.figure(figsize=(10, 10))
    plt.xlim(-200, 200)
    plt.ylim(-200, 200)
    # tmp = np.ma.masked_where(nc['ref_sum'][:] == 0, nc['ref_sum'][:])
    plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, pow2db(nc['ref_linear_sum'][:] / nc.num_scans),
                   vmin=0, cmap='nipy_spectral')
    plt.colorbar()
    contours = plt.contour(x_m * 1e-3, y_m * 1e-3, pow2db(nc['ref_linear_sum'][:] / nc.num_scans),
                           levels=[np.max(pow2db(nc['ref_linear_sum'][:] / nc.num_scans))/7.5])
    plt.title("Average Reflectivity " + nc_file.split('\\')[-1].replace('.nc', ''))

    # Significant bat areas to save
    sbis = []

    # Loop through the contours
    for path in contours.collections[0].get_paths():
        try:
            poly = Polygon(path.vertices)
            if MIN_AREA < poly.area < MAX_AREA:
                # Get the xy coord of the verticies
                x, y = poly.exterior.xy

                # plot these points
                plt.plot(x, y, color='r')

                # Get the representative point of the polygon
                rep_x = poly.representative_point().x
                rep_y = poly.representative_point().y

                # get the area of the point
                area = poly.area

                # Create dictionary with required information and append it to the list
                poly_dict = {'x': x,
                             'y': y,
                             'area': area,
                             'rep_x': rep_x,
                             'rep_y': rep_y
                             }

                sbis.append(poly_dict)
        except Exception as e:
            print("Exception Encountered! Skipping the contour")
            print(e)
            pass

    # Get the cave info for this radar
    caves = read_caves(cave_csv, nc.site)

    # Convert cave lat/lon to x/y coord system
    radar_lat = nc.radar_lat
    radar_lon = nc.radar_lon
    x_radar, y_radar, _, _ = utm.from_latlon(radar_lat, radar_lon)

    for cave in caves:
        # Convert lat-lons to utm for
        x_bat, y_bat, _, _ = utm.from_latlon(float(cave['lat']), float(cave['lon']))

        # Calc relative x and y of roost
        x_rel_m = (x_bat - x_radar)
        y_rel_m = (y_bat - y_radar)

        cave['x'] = x_rel_m/1e3
        cave['y'] = y_rel_m/1e3
        plt.scatter(cave['x'], cave['y'], c='y')

        # Get rid of misc crap
        del x_bat, y_bat, _

    # Find the closest cave theo the sbi
    for sbi in sbis:
        closest = None
        for i, cave in enumerate(caves):
            dist = distance(cave['x'], sbi['rep_x'], cave['y'], sbi['rep_y'])
            # print(closest, dist)
            if closest is None and dist < MAX_DISTANCE:
                closest = (dist, i)
            elif MAX_DISTANCE > dist < closest[0]:
                closest = (dist, i)
            plt.scatter(cave['x'], cave['y'], c='y')
        if closest is not None:
            sbi['cave'] = caves[closest[1]]['name']
            sbi['cave_x'] = caves[closest[1]]['x']
            sbi['cave_y'] = caves[closest[1]]['y']
            plt.plot([sbi['cave_x'], sbi['rep_x']], [sbi['cave_y'], sbi['rep_y']], color='g')

        plt.scatter(sbi['rep_x'], sbi['rep_y'])

    plt.savefig(os.path.join('data\\images8', nc_file.split('\\')[-1].replace('.nc', '.png')))
    plt.close()
    nc.close()


if __name__ == '__main__':
    files = glob('data/netcdf_daily_average/K*.nc')
    for f in files:
        print(f)
        get_sba(f, 'cave_locations.csv')
        # try:
        #     get_sba(f, 'cave_locations.csv')
        # except Exception as e:
        #     print(e)
        #     continue
