import matplotlib.pyplot as plt
import numpy as np


def db2pow(x):
    return 10.**(x/10.)


def pow2db(x):
    return 10. * np.log10(x)


def save_data_as_image(x, y, data, image_name):
    fig = plt.figure(figsize=(10, 10), dpi=100)
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    plt.pcolormesh(x, y, data, cmap='nipy_spectral')
    plt.axis('off')
    fig.get_axes()[0].get_xaxis().set_visible(False)
    fig.get_axes()[0].get_yaxis().set_visible(False)
    plt.savefig(image_name, bbox_inches='tight', pad_inches=0)


# Lat lon of Frio cave
lat_bat = 29.434507
lon_bat = -99.684534