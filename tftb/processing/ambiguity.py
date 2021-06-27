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
from scipy.signal import hilbert
from tftb.utils import nextpow2


def _find(condt):
    res, = np.nonzero(np.ravel(condt))
    return res


def wide_band(signal, fmin=None, fmax=None, N=None):
    if 1 in signal.shape:
        signal = signal.ravel()
    elif signal.ndim != 1:
        raise ValueError("The input signal should be one dimensional.")
    s_ana = hilbert(np.real(signal))
    nx = signal.shape[0]
    m = int(np.round(nx / 2.0))
    t = np.arange(nx) - m
    tmin = 0
    tmax = nx - 1
    T = tmax - tmin

    # determine default values for fmin, fmax
    if (fmin is None) or (fmax is None):
        STF = np.fft.fftshift(s_ana)
        sp = np.abs(STF[:m]) ** 2
        maxsp = np.amax(sp)
        f = np.linspace(0, 0.5, m + 1)
        f = f[:m]
        indmin = _find(sp > maxsp / 100.0).min()
        indmax = _find(sp > maxsp / 100.0).max()
        if fmin is None:
            fmin = max([0.01, 0.05 * np.fix(f[indmin] / 0.05)])
        if fmax is None:
            fmax = 0.05 * np.ceil(f[indmax] / 0.05)
    B = fmax - fmin
    R = B / ((fmin + fmax) / 2.0)
    nq = np.ceil((B * T * (1 + 2.0 / R) * np.log((1 + R / 2.0) / (1 - R / 2.0))) / 2.0)  # NOQA
    nmin = nq - (nq % 2)
    if N is None:
        N = int(2 ** (nextpow2(nmin)))

    # geometric sampling for the analyzed spectrum
    k = np.arange(1, N + 1)
    q = (fmax / fmin) ** (1.0 / (N - 1))
    geo_f = fmin * (np.exp((k - 1) * np.log(q)))
    tfmatx = -2j * np.dot(t.reshape(-1, 1), geo_f.reshape(1, -1)) * np.pi
    tfmatx = np.exp(tfmatx)
    S = np.dot(s_ana.reshape(1, -1), tfmatx)
    S = np.tile(S, (nx, 1))
    Sb = S * tfmatx

    tau = t
    S = np.c_[S, np.zeros((nx, N))].T
    Sb = np.c_[Sb, np.zeros((nx, N))].T

    # mellin transform computation of the analyzed signal
    p = np.arange(2 * N)
    coef = np.exp(p / 2.0 * np.log(q))
    mellinS = np.fft.fftshift(np.fft.ifft(S[:, 0] * coef))
    mellinS = np.tile(mellinS, (nx, 1)).T

    mellinSb = np.zeros((2 * N, nx), dtype=complex)
    for i in range(nx):
        mellinSb[:, i] = np.fft.fftshift(np.fft.ifft(Sb[:, i] * coef))

    k = np.arange(1, 2 * N + 1)
    scale = np.logspace(np.log10(fmin / fmax), np.log10(fmax / fmin), N)
    theta = np.log(scale)
    mellinSSb = mellinS * np.conj(mellinSb)

    waf = np.fft.ifft(mellinSSb, N, axis=0)
    no2 = int((N + N % 2) / 2.0)
    waf = np.r_[waf[no2:(N + 1), :], waf[:no2, :]]

    # normalization
    s = np.real(s_ana)
    SP = np.fft.fft(hilbert(s))
    indmin = int(1 + np.round(fmin * (nx - 2)))
    indmax = int(1 + np.round(fmax * (nx - 2)))
    sp_ana = SP[(indmin - 1):indmax]
    waf *= (np.linalg.norm(sp_ana) ** 2) / waf[no2 - 1, m - 1] / N

    return waf, tau, theta


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
        n_fbins = signal.shape[0]

    naf = np.zeros((n_fbins, taucol), dtype=complex)
    for icol in range(taucol):
        taui = int(lag[icol])
        t = np.arange(abs(taui), n - abs(taui)).astype(int)
        naf[t, icol] = signal[t + taui] * np.conj(signal[t - taui])
    naf = np.fft.fft(naf, axis=0)

    _ix1 = np.arange((n_fbins + (n_fbins % 2)) // 2, n_fbins)
    _ix2 = np.arange((n_fbins + (n_fbins % 2)) // 2)

    _xi1 = -(n_fbins - (n_fbins % 2)) // 2
    _xi2 = ((n_fbins + (n_fbins % 2)) // 2 - 1)
    xi = np.arange(_xi1, _xi2 + 1, dtype=float) / n_fbins
    naf = naf[np.hstack((_ix1, _ix2)), :]
    return naf, lag, xi


if __name__ == '__main__':
    from tftb.generators.misc import altes
    sig = altes(128, 0.1, 0.45)
    waf, tau, theta = wide_band(sig)
    from matplotlib.pyplot import contour, show
    contour(tau, theta, np.abs(waf) ** 2)
    show()
