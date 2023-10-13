__all__ = [
    "abs2",
    "butt_low",
    "sav_smooth",
    "extract_flucts",
    "fluct_rms",
    "calc_psd",
    "calc_xspec",
]
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt, savgol_filter
from scipy.optimize import curve_fit


def abs2(x):
    """
    Calculates the absolute square of a complex number. Works on ND arrays.

    Args:
        x (_type_): complex number or ND array of complex numbers

    Returns:
        _type_: absolute square of x
    """
    a = x.real**2 + x.imag**2
    return a


def butt_low(sig, cut_freq, samp_freq, order=4):
    """
    Wrapper for scipy.signal.butter to filter a signal with a low pass filter.
    Default order is 4. Works on 1D and multi-D arrays.
    For multi-D arrays, the last dimension is assumed to be the time dimension.

    Args:
        sig (_type_): ND array of signal to be filtered
        cut_freq (_type_): cutoff frequency
        samp_freq (_type_): sampling frequency
        order (int, optional): Butterworth filter order. Defaults to 4.

    Returns:
        _type_: filtered signal array
    """
    sos = butter(order, cut_freq, "lp", fs=samp_freq, output="sos")
    filtered_sig = sosfiltfilt(sos, sig)
    return filtered_sig


def sav_smooth(a, b=10, axis=-1):
    """
    Wrapper for scipy.signal.savgol_filter to smooth a signal. Default window size is 1/10th of the signal length.
    Works on 1D and multi-D arrays (maybe? untested).

    Args:
        a (_type_): Array to be smoothed
        b (int, optional): Number of equal sized non overlapping windows. Defaults to 10.

    Returns:
        _type_: smoothed array
    """
    wind = a.shape[-1] // (b)
    if wind % 2 == 0:
        wind += 1
    c = savgol_filter(a, wind, 1, axis=axis)
    return c


def extract_flucts(sigs, smooth_sigs):
    """
    Extracts fluctuations from a signal by subtracting a smoothed version of the signal from the original signal. Works on 1D and multi-D arrays.
    """
    return sigs - smooth_sigs


def fluct_rms(flucts):
    """
    Calculates the root mean square of fluctuations. Works on 1D and multi-D arrays.

    Args:
        flucts (_type_): Array of fluctuations (signal - smoothed signal)

    Returns:
        _type_: root mean square of fluctuations
    """
    return np.sqrt(np.mean(flucts**2, axis=-1))


def calc_psd(a):
    """
    Calculates the power spectral density of a signal. Works on 1D and multi-D arrays. (verify if a real fft should be used instead of a complex fft)

    Args:
        a (_type_): _description_

    Returns:
        _type_: _description_
    """
    afft = np.fft.rfft(a, axis=-1)
    return abs2(afft)


def calc_xspec(a, b):
    """
    Calculates the cross spectral density of two signals. Works on 1D and multi-D arrays. (check what type of input is needed for this function)

    Args:
        a (_type_): _description_
        b (_type_): _description_

    Returns:
        _type_: _description_
    """
    return a * np.conj(b)


# def detrend_2d(sigs):
#     filt = butt_low(sigs, 1e5, fs, order=4)[:, :t_end]
#     avg = np.mean(filt.reshape(nx * ny, nshots, t_end), axis=1)
#     pt = np.mean(avg[:, t_end - 10000 : t_end], axis=-1)
#     smooth = sav_smooth(avg, t_end)
#     flucts = extract_flucts(
#         np.moveaxis(filt.reshape(nx * ny, nshots, t_end), 1, 0), smooth
#     )
#     rms = fluct_rms(flucts)
#     return filt, pt, flucts, rms


# def remove_offset(a, b):
#     offsets = np.min(a)
#     c = a - offsets
#     return c


# def remove_offset2(a, b):
#     offsets = np.mean(a[:, :, b:-1], axis=-1)
#     c = np.moveaxis((np.moveaxis(a, -1, 0) - offsets), 0, -1)
#     return c


# def isat_to_density(isat, temp, p_area):
#     dens = (isat * 1e-3 * 1e-6) / (np.sqrt(k * evK * temp / M) * 0.61 * p_area * echar)
#     return dens


# def lapd_xy_grid(nx, dx, x0, ny, dy, y0):
#     """
#     Generates a 2D grid of x and y coordinates using the LAPD coordinate system for Bz out of the page.

#     Args:
#         nx (_type_): _description_
#         dx (_type_): _description_
#         x0 (_type_): _description_
#         ny (_type_): _description_
#         dy (_type_): _description_
#         y0 (_type_): _description_

#     Returns:
#         _type_: _description_
#     """
#     x = np.linspace(x0 - (nx // 2) * dx, x0 + (nx // 2) * dx, nx, endpoint=True)
#     y = np.linspace(y0 + (ny // 2) * dy, y0 - (ny // 2) * dy, ny, endpoint=True)
#     return x, y, np.meshgrid(x, y)
