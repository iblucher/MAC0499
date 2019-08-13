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
    if filename != None:
        title = filename.split('/')[-1].split('.')[0]
    else:
        title = filename
    
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.plot(wavelength, flux, '-')
    plt.xlabel(r"Wavelength")
    plt.ylabel(r"Flux")
    plt.title(title)
    plt.show()


def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    print(idx - 1, idx)
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def fast_dtw_on_stellar_spectra(tel, obs, dist=euclidean):
    distance, path = fastdtw(tel, obs, dist=dist)
    
    print('Distance: {}'.format(distance))
    
    warp_idxs = all(tel == obs for tel, obs in path)
    print('Are all warp path indices aligned? {}'.format(warp_idxs))
    
    fast_path = list(zip(*path))
    plt.plot(fast_path[0], fast_path[1], 'k')
    plt.show()
    
    return fast_path


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