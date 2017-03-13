import csv
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import os
import utm
from shapely.geometry import Polygon
from glob import glob
from utils import greens_theorem, read_caves


def distance(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2 - y1)**2)


def get_sba(nc_file, cave_csv):
    nc = netCDF4.Dataset(nc_file)

    # Get the grid figured out
    range_m, az_rad = np.meshgrid(nc['range'][:], nc['azimuth'][:])
    elev = np.deg2rad(nc.elevation)

    x_m = range_m * np.cos(elev) * np.sin(az_rad)
    y_m = range_m * np.cos(elev) * np.cos(az_rad)

    plt.figure(figsize=(10, 10))
    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, nc['phi_dp_weighted_sum'][:] / nc.num_scans,
                   vmin=0, vmax=360, cmap='nipy_spectral')
    contours = plt.contour(x_m * 1e-3, y_m * 1e-3, nc['phi_dp_weighted_sum'][:] / nc.num_scans,
                           levels=[360.])
    plt.title("PhiDP")

    # Significant bat areas to save
    sbis = []

    # Loop through the contours
    for path in contours.collections[0].get_paths():
        verts = path.vertices

        if greens_theorem(verts) > 100:
            try:
                # Convert the contour to a polygon
                poly = Polygon(verts)

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
            except Exception:
                print("Exception Encountered! Skipping the contour")
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

        # Get rid of misc crap
        del x_bat, y_bat, _

    # Find the closest cave to the sbi
    for sbi in sbis:
        closest = None
        for i, cave in enumerate(caves):
            dist = distance(cave['x'], sbi['rep_x'], cave['y'], sbi['rep_y'])
            print(dist)
            if closest is None:
                closest = (dist, i)
            elif dist < closest[0]:
                closest = (dist, i)
            print(cave['x'], cave['y'])
            plt.scatter(cave['x'], cave['y'], c='y')


        print(closest)

        sbi['cave'] = caves[closest[1]]['name']
        sbi['cave_x'] = caves[closest[1]]['x']
        sbi['cave_y'] = caves[closest[1]]['y']

        plt.scatter(sbi['rep_x'], sbi['rep_y'])

    plt.show()


if __name__ == '__main__':
    get_sba('/Users/tbupper90/data/bci_study/testing/netcdf/average_20160617.nc',
            '/Users/tbupper90/Google Drive/Research/BCI_study/cave_locations.csv')
