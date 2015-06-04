#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Ambiguity functions.
"""

import numpy as np


def narrow_band(signal, lag=None, n_fbins=None):
    """Narrow band ambiguity function.

    :param signal: Signal to be analyzed.
    :param lag: vector of lag values.
    :param n_fbins: number of frequency bins
    :type signal: array-like
    :type lag: array-like
    :type n_fbins: int
    :return: Doppler lag representation
    :rtype: array-like
    """

    n = signal.shape[0]
    if lag is None:
        if n % 2 == 0:
            tau_start, tau_end = -n / 2 + 1, n / 2
        else:
            tau_start, tau_end = -(n - 1) / 2, (n + 1) / 2
        lag = np.arange(tau_start, tau_end)
    taucol = lag.shape[0]

    if n_fbins is None:
        n_fbins = n

    naf = np.zeros((n_fbins, taucol), dtype=complex)
    for icol in xrange(taucol):
        taui = lag[icol]
        t = np.arange(abs(taui), n - abs(taui))
        naf[t, icol] = signal[t + taui] * np.conj(signal[t - taui])
    naf = np.fft.fft(naf, axis=0)

    _ix1 = np.arange((n_fbins + (n_fbins % 2)) / 2, n_fbins)
    _ix2 = np.arange((n_fbins + (n_fbins % 2)) / 2)

    _xi1 = -(n_fbins - (n_fbins % 2)) / 2
    _xi2 = ((n_fbins + (n_fbins % 2)) / 2 - 1)
    xi = np.arange(_xi1, _xi2 + 1, dtype=float) / n_fbins
    naf = naf[np.hstack((_ix1, _ix2)), :]
    return naf, lag, xi

if __name__ == '__main__':
    from tftb.generators.api import anabpsk
    sig = anabpsk(256, 8)[0]
    naf, x, y = narrow_band(sig)
    import matplotlib.pyplot as plt
    plt.contour(2 * x, y, np.abs(naf) ** 2, 16)
    plt.grid(True)
    plt.show()
