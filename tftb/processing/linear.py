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
        try:
            tfr[indices.astype(int), icol] = signal[ti + tau - 1] * np.conj(window[lh + tau])
        except IndexError:
            from IPython.core.debugger import Tracer
            Tracer()()
    tfr = np.fft.fft(tfr, axis=0)

    if n_fbins % 2 == 0:
        freqs = np.hstack((np.arange(n_fbins / 2), np.arange(-n_fbins / 2, 0)))
    else:
        freqs = np.hstack((np.arange((n_fbins - 1) / 2),
                           np.arange(-(n_fbins - 1) / 2, 0)))
    return tfr, time_instants, freqs

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from tftb.generators.api import fmconst
    sig = np.hstack((fmconst(128, 0.2)[0], fmconst(128, 0.4)[0]))
    tfr, t, f = stft(sig)
    plt.imshow(np.abs(tfr))
    plt.show()
