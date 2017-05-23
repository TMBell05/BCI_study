import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import csv
from scipy.interpolate import RegularGridInterpolator


def interpolate(az, rng, data, new_az, new_rng):

    # TODO - Do this better.....
    # Quick fix for repeated azimuths
    tmp, cts = np.unique(az, return_counts=True)
    repeated_values = tmp[cts > 1]
    for value in repeated_values:
        ind = np.where(az == value)
        az[ind[0][1:]] += .0001


    # Make the interpolate function
    interp_funct = RegularGridInterpolator((az, rng), data)
    interp_funct.bounds_error = False

    # Make the points
    pts = np.asarray([[a, r] for a, r in zip(new_az.flatten(), new_rng.flatten())])

    # do the interpolation
    new_data = interp_funct(pts)
    return new_data.reshape(new_az.shape)


def db2pow(x):
    return 10.**(x/10.)


def distance(x1, x2, y1, y2):
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)


def pow2db(x):
    return 10. * np.log10(x)


def read_caves(csv_file, site):
    caves = []

    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['radar'] == site:
                caves.append(row)

    return caves


def save_data_as_tiff(x, y, data, image_name):
    fig = plt.figure(figsize=(10, 10), dpi=150)
    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.pcolormesh(x, y, data, cmap='gray', vmin=0, vmax=360)
    plt.axis('off')
    fig.get_axes()[0].get_xaxis().set_visible(False)
    fig.get_axes()[0].get_yaxis().set_visible(False)
    plt.savefig(image_name, bbox_inches='tight', pad_inches=0)
    plt.close()


def greens_theorem(vs):
    a = 0
    x0,y0 = vs[0]
    for [x1,y1] in vs[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a


def plot_netcdf(netcdf_file, img_name=None):
    nc = netCDF4.Dataset(netcdf_file)

    range_m, az_rad = np.meshgrid(nc['range'][:], nc['azimuth'][:])
    # az_rad = np.deg2rad(az_deg)
    elev = np.deg2rad(nc.elevation)

    x_m = range_m * np.cos(elev) * np.sin(az_rad)
    y_m = range_m * np.cos(elev) * np.cos(az_rad)

    plt.figure(figsize=(16, 8))
    plt.subplot(1, 2, 1)
    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, nc['phi_dp_sum'][:]/nc.num_scans, vmin=0, vmax=360, cmap='nipy_spectral')
    plt.colorbar()
    plt.title("PhiDP")

    plt.subplot(1, 2, 2)
    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, nc['phi_dp_weighted_sum'][:]/nc.num_scans, vmin=0, vmax=360, cmap='nipy_spectral')
    plt.colorbar()
    plt.title("PhiDP Weighted")

    # plt.subplot(2, 2, 3)
    # plt.xlim(-150, 150)
    # plt.ylim(-150, 150)
    # plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, nc['ref_linear_sum'][:] / nc.num_scans)
    # plt.colorbar()
    # plt.title("Reflectivity")
    #
    # plt.subplot(2, 2, 4)
    # plt.xlim(-150, 150)
    # plt.ylim(-150, 150)
    # plt.pcolormesh(x_m * 1e-3, y_m * 1e-3, pow2db(nc['eta_linear_sum'][:] / nc.num_scans))
    # plt.colorbar()
    # plt.title("Eta")

    # plt.suptitle('KDFX %Y%m%d-%H%M%S Elev: {elev}'.format(elev=np.rad2deg(elev)))
    # plt.close()
    if img_name is not None:
        plt.savefig(img_name)
        save_data_as_tiff(x_m/1e3, y_m/1e3, nc['phi_dp_weighted_sum'][:]/nc.num_scans, 'test.png')

    nc.close()
    plt.close()
    # plt.show()

# Lat lon of Frio cave
lat_bat = 29.434507
lon_bat = -99.684534