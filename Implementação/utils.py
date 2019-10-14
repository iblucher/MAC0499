import math
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

import astropy.io.fits as pf
from fastdtw import fastdtw


euclidean = euclidean_norm = lambda x, y: np.abs(x - y)


def get_spectrum_data_and_wavelength(file):
    spec_data = pf.getdata(file)
    spec_header = pf.getheader(file)
    
    wl_i = spec_header['CRVAL1'] # Coordinate value of position in degrees, specified in CRPIX.
    wl_step = spec_header['CDELT1'] # Increment-per-pixel of axis n, in degrees
    wavelength = np.arange(spec_data.size)*wl_step+wl_i
    
    return spec_data, wavelength


def plot_spectrum(flux, wavelength, filename=None):
    title = filename
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.plot(wavelength, flux, '-')
    plt.xlabel(r"Wavelength (nm)")
    plt.ylabel(r"Flux")
    plt.title(title)
    plt.show()


def find_nearest(array,value):
    idx = np.searchsorted(array, value)
    return(idx - 1, idx)

def fast_dtw_on_stellar_spectra(tel, obs, dist=euclidean, xlabel=None, ylabel=None):
    distance, path = fastdtw(tel, obs, dist=dist)
    
    print('Distance: {}'.format(distance))
    
    warp_idxs = all(tel == obs for tel, obs in path)
    print('Are all warp path indices aligned? {}'.format(warp_idxs))
    
    fast_path = list(zip(*path))
    plt.plot(fast_path[0], fast_path[1], 'k')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    
    return fast_path, path


def align_sequence_dtw_path(path, original_seq):
    zipped_path = zip(path[0], path[1])
    
    path_dict = defaultdict(list)
    for i, j in zipped_path:
        path_dict[j].append(i)
        
    aligned_seq = []    
    for i, j in path_dict.items():
        same_path_idx = [original_seq[t] for t in j]
        aligned_seq.insert(i, np.mean(same_path_idx))
        
    return aligned_seq