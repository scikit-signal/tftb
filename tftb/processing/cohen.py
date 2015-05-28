#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Bilinear Time-Frequency Processing in the Cohen’s Class.
"""

import warnings
import numpy as np

from tftb.utils import nextpow2


def spectrogram(signal, time_instants=None, n_fbins=None, window=None):
    """Compute the spectrogram of a signal.

    :param signal: Signal to be analyzed.
    :param time_instants: timestamps of the signal.
    :param n_fbins: number of frequency bins
    :param window: analysis window
    :type signal: array-like
    :type time_instants: array-like
    :type n_fbins: int
    :type window: array-like
    :return: time frequency representation
    :rtype: array-like
    """
    if time_instants is None:
        time_instants = np.arange(signal.shape[0])

    if n_fbins is None:
        n_fbins = signal.shape[0]

    if window is None:
        hlength = np.floor(signal.shape[0] / 4.0)
        hlength += 1 - np.remainder(hlength, 2)
        from scipy.signal import hamming
        window = hamming(hlength)
    else:
        hlength = window.shape[0]
        if hlength % 2 == 0:
            raise ValueError("Smoothing window should have an odd length.")

    if 2 ** nextpow2(n_fbins) != n_fbins:
        msg = "For faster computation, the frequency bins should be a power of 2."
        warnings.warn(msg, UserWarning)

    lh = (window.shape[0] - 1) / 2
    tfr = np.zeros((n_fbins, time_instants.shape[0]), dtype=complex)
    for icol in xrange(time_instants.shape[0]):
        ti = time_instants[icol]
        start = min([np.round(n_fbins / 2.0) - 1, lh, ti - 1])
        end = min([np.round(n_fbins / 2.0) - 1, lh, signal.shape[0] - ti - 1])
        tau = np.arange(-start, end)
        indices = np.remainder(n_fbins + tau, n_fbins).astype(int)
        tfr[indices, icol] = signal[ti + tau] * np.conj(window[lh + tau]) / np.linalg.norm(window[lh + tau])

    if n_fbins % 2 == 0:
        left = np.arange(n_fbins / 2, dtype=float) / n_fbins
        right = np.arange(-n_fbins / 2, 0, dtype=float) / n_fbins
    else:
        left = np.arange((n_fbins - 1) / 2, dtype=float) / n_fbins
        right = np.arange(-(n_fbins - 1) / 2, 0, dtype=float) / n_fbins
    frequencies = np.hstack((left, right))

    return np.abs(np.fft.fft(tfr, axis=0)) ** 2, time_instants, frequencies


def pseudo_wigner_ville(signal, time_samples=None, freq_bins=None, window=None):
    """Compute the Pseudo Wigner Ville time frequency distribution.

    :param signal: Signal to analyze.
    :param time_samples: time instants
    :param freq_bins: Number of frequency bins.
    :param window: Frequency smoothing window (Default: Hamming of length
    signal-length / 4)
    :type signal: array-like.
    :type time_samples: array-like
    :type freq_bins: int
    :type window: array-like
    :return: Pseudo Wigner Ville time frequency distribution.
    :rtype: array-like
    """
    if time_samples is None:
        time_samples = np.arange(signal.shape[0])

    if freq_bins is None:
        freq_bins = signal.shape[0]

    if 2 ** nextpow2(freq_bins) != freq_bins:
        msg = "For faster computation, the frequency bins should be a power of 2."
        warnings.warn(msg, UserWarning)

    if window is None:
        winlength = np.floor(signal.shape[0] / 4.0)
        winlength = winlength + 1 - np.remainder(winlength, 2)
        from scipy.signal import hamming
        window = hamming(int(winlength))
    elif window.shape[0] % 2 == 0:
        raise ValueError('The smoothing window must have an odd length.')

    tfr = np.zeros((freq_bins, time_samples.shape[0]), dtype=complex)
    lh = (window.shape[0] - 1) / 2
    for icol in xrange(time_samples.shape[0]):
        ti = time_samples[icol]
        taumaxvals = (ti, signal.shape[0] - ti - 1,
                      np.round(freq_bins / 2.0), lh)
        taumax = np.min(taumaxvals)
        tau = np.arange(-taumax, taumax + 1).astype(int)
        indices = np.remainder(freq_bins + tau, freq_bins).astype(int)
        tfr[indices, icol] = window[lh + tau] * signal[ti + tau] * \
            np.conj(signal[ti - tau])
        tau = np.round(freq_bins / 2.0)
        if (ti <= signal.shape[0] - tau) and (ti >= tau + 1) and (tau <= lh):
            tfr[int(tau), icol] = 0.5 * (window[lh + tau] * signal[ti + tau, 0] *
                    np.conj(signal[ti - tau, 0]) + window[lh - tau] *
                    signal[ti - tau, 0] * np.conj(signal[ti + tau, 0]))

    tfr = np.fft.fft(tfr, axis=0)
    return np.real(tfr)


def wigner_ville(signal, time_samples=None, freq_bins=None):
    """Compute the Wigner Ville time frequency distribution.

    :param signal: Signal to analyze.
    :param time_samples: time instants
    :param freq_bins: Number of frequency bins.
    :type signal: array-like.
    :type time_samples: array-like
    :type freq_bins: int
    :return: Wigner Ville time frequency distribution.
    :rtype: array-like
    """
    if time_samples is None:
        time_samples = np.arange(signal.shape[0])

    if freq_bins is None:
        freq_bins = signal.shape[0]

    if 2 ** nextpow2(freq_bins) != freq_bins:
        msg = "For faster computation, the frequency bins should be a power of 2."
        warnings.warn(msg, UserWarning)

    tfr = np.zeros((freq_bins, time_samples.shape[0]), dtype=complex)

    for icol in xrange(time_samples.shape[0]):
        ti = time_samples[icol]
        taumax = np.min((ti, signal.shape[0] - ti - 1,
                         np.round(freq_bins / 2.0) - 1))
        tau = np.arange(-taumax, taumax + 1).astype(int)
        indices = np.remainder(freq_bins + tau, freq_bins).astype(int)
        tfr[indices, icol] = signal[ti + tau] * np.conj(signal[ti - tau])
        tau = np.round(freq_bins / 2.0)
        if (ti <= signal.shape[0] - tau) and (ti >= tau + 1):
            tfr[tau, icol] = 0.5 * (signal[ti + tau, 0] * np.conj(signal[ti - tau, 0])) + \
                                   (signal[ti - tau, 0] * np.conj(signal[ti + tau, 0]))

    tfr = np.fft.fft(tfr, axis=0)
    tfr = np.real(tfr)
    return tfr


if __name__ == '__main__':
    from tftb.generators.api import fmlin
    signal, _ = fmlin(128, 0.1, 0.4)
    from scipy.signal import kaiser
    window = kaiser(17, 3 * np.pi)
    tfr, time, freq = spectrogram(signal, n_fbins=64, window=window)
    import matplotlib.pyplot as plt
    plt.imshow(tfr)
    plt.show()
