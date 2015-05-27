#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Linear Time Frequency Processing.
"""

import numpy as np
from tftb.utils import nearest_odd, divider, modulo, izak


def gabor(signal, n_coeff=None, q_oversample=None, window=None):
    """Compute the Gabor representation of a signal.

    :param signal: Singal to be analyzed.
    :param n_coeff: number of Gabor coefficients in time.
    :param q_oversample: Degree of oversampling
    :param window: Synthesis window
    :type signal: array-like
    :type n_coeff: integer
    :type q_oversample: int
    :type window: array-like
    :return: Tuple of Gabor coefficients, biorthogonal window associated with the synthesis window.
    :rtype: tuple
    """
    if n_coeff is None:
        n_coeff = divider(signal.shape[0])
    if q_oversample is None:
        q_oversample = divider(n_coeff)
    if window is None:
        window = np.exp(np.log(0.005) * np.linspace(-1, 1, nearest_odd(n_coeff)) ** 2)
        window = window / np.linalg.norm(window)
    m = q_oversample * signal.shape[0] / float(n_coeff)
    mb = signal.shape[0] / float(n_coeff)
    nb = signal.shape[0] / float(m)

    # Zak transform?
    nh = window.shape[0]
    if nh % 2 == 0:
        raise ValueError("The window function should have an odd length.")
    alpha = np.round((2 * signal.shape[0] / float(n_coeff) - 1 - nh) / (2 *
        q_oversample))
    hn1 = np.zeros((signal.shape[0],))
    start = np.round(((signal.shape[0] - (nh - 1))) / 2) - alpha
    end = np.round((signal.shape[0] + nh - 1) / 2) - alpha
    hn1[np.arange(start - 1, end).astype(int)] = window

    msig = hn1.reshape(nb, m, order='F')
    dzth = np.fft.fft(msig.T, axis=0) / np.sqrt(m)
    mzh = np.zeros((m, mb))
    x = np.arange(1, m + 1, dtype=float)
    for l in xrange(q_oversample):
        mod = modulo(x - l * m / q_oversample, m).astype(int)
        mzh += np.abs(dzth[mod - 1, :]) ** 2

    mzh[mzh < np.spacing(1)] = 1

    # Za transform of biorthogonal dual frame window gam
    dztgam = dzth / mzh
    gam = np.real(izak(dztgam)) / signal.shape[0]

    # Computation of Gabor coefficient of dual frame window.
    dgrn1 = np.zeros((signal.shape[0], n_coeff), dtype=complex)
    k = np.arange(1, signal.shape[0] + 1)
    for n in xrange(n_coeff):
        index = modulo(k - n * m / q_oversample, signal.shape[0]).astype(int) - 1
        dgrn1[:, n] = np.fft.fft(signal * np.fft.fftshift(gam[index]), axis=0)
    dgr = dgrn1[np.arange(signal.shape[0], step=nb).astype(int), :]
    tfr = np.abs(dgr) ** 2
    return tfr, dgr, gam


def stft(signal, time_instants=None, n_fbins=None, window=None):
    """Compute the short time Fourier transform.

    :param signal: signal to be analyzed.
    :param time_instants: timestamps of the signal.
    :param n_fbins: Number of frequency bins.
    :param window: Window function for frequency smoothing.
    :type signal: array-like
    :type time_instants: array-like
    :type n_fbins: int
    :type window: array-like
    :return: tuple of (tfr, timestamps, normalized frequency vector)
    :rtype: tuple
    """
    signal = signal.ravel()
    if time_instants is None:
        time_instants = np.arange(1, signal.shape[0])
    if n_fbins is None:
        n_fbins = signal.shape[0]
    if window is None:
        hlength = np.floor(n_fbins / 4.0)
        hlength = hlength + 1 - np.remainder(hlength, 2)
        from scipy.signal import hamming
        window = hamming(int(hlength))

    lh = (window.shape[0] - 1) / 2
    window = window / np.linalg.norm(window)

    tfr = np.zeros((n_fbins, time_instants.shape[0]), dtype=complex)
    for icol in xrange(tfr.shape[1]):
        ti = time_instants[icol]
        start = -np.min([np.round(n_fbins / 2.0) - 1, lh, ti - 1])
        end = np.min([np.round(n_fbins / 2.0) - 1, lh, signal.shape[0] - ti])
        tau = np.arange(start, end + 1).astype(int)
        indices = np.remainder(n_fbins + tau, n_fbins)
        tfr[indices.astype(int), icol] = signal[ti + tau - 1] * np.conj(window[lh + tau])
    tfr = np.fft.fft(tfr, axis=0)

    if n_fbins % 2 == 0:
        freqs = np.hstack((np.arange(n_fbins / 2), np.arange(-n_fbins / 2, 0)))
    else:
        freqs = np.hstack((np.arange((n_fbins - 1) / 2),
                           np.arange(-(n_fbins - 1) / 2, 0)))
    return tfr, time_instants, freqs

if __name__ == '__main__':
    from tftb.generators.api import fmlin
    sig, _ = fmlin(128)
    a, b, c = gabor(sig, 64, 32)
    from matplotlib.pyplot import imshow, show
    imshow(a)
    show()
