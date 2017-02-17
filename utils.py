import matplotlib.pyplot as plt
import numpy as np
import netCDF4


def db2pow(x):
    return 10.**(x/10.)


def pow2db(x):
    return 10. * np.log10(x)


def save_data_as_tiff(x, y, data, image_name):
    fig = plt.figure(figsize=(10, 10), dpi=150)
    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.pcolormesh(x, y, data, cmap='gray', vmin=0, vmax=360)
    plt.axis('off')
    fig.get_axes()[0].get_xaxis().set_visible(False)
    fig.get_axes()[0].get_yaxis().set_visible(False)
    plt.savefig(image_name, bbox_inches='tight', pad_inches=0)


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
        save_data_as_tiff(x_m/1e3, y_m/1e3, nc['phi_dp_weighted_sum'][:]/nc.num_scans, 'test.png')
        plt.savefig(img_name)
    plt.close()
    # plt.show()

# Lat lon of Frio cave
lat_bat = 29.434507
lon_bat = -99.684534