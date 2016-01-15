#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Bilinear Time-Frequency Processing in the Affine Class.
"""

import numpy as np
from matplotlib.mlab import find
from scipy.signal import hilbert
from scipy.optimize import brenth, newton
from scipy.interpolate import splrep, splev
from tftb.generators import mexhat, scale
from tftb.processing.utils import integrate_2d
from tftb.utils import nextpow2


def lambdak(u, k):
    method_lookup = {'bertrand': 0, 'unterberger': -1, 'd_flandrin': 0.5,
                     'aspwv': 2}
    k = method_lookup[k]
    if k not in (0, 1.0):
        y = (k * (np.exp(-u) - 1.0) / (np.exp(-k * u) - 1.0))
        y = y ** (1.0 / (k - 1.0))
    elif k == 1.0:
        y = np.exp(1.0 + u * np.exp(-u) / (np.exp(-u) - 1.0))
    elif k == 0:
        if u == 0:
            return 1.0
        y = -u / (np.exp(-u) - 1.0)
    return y


def smoothed_pseudo_wigner(signal, timestamps=None, K='bertrand', nh0=None,
        ng0=0, fmin=None, fmax=None, n_voices=None):
    """smoothed_pseudo_wigner

    :param signal:
    :param timestamps:
    :param K:
    :param nh0:
    :param ng0:
    :param fmin:
    :param fmax:
    :param n_voices:
    :type signal:
    :type timestamps:
    :type K:
    :type nh0:
    :type ng0:
    :type fmin:
    :type fmax:
    :type n_voices:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(signal.shape[0])
    if nh0 is None:
        nh0 = np.sqrt(signal.shape[0])

    tcol = timestamps.shape[0]
    mt = signal.shape[0]

    x1 = x2 = signal.copy()
    s1 = np.real(x1)
    s2 = np.real(x2)
    m = (mt + np.remainder(mt, 2.0)) / 2.0

    if (fmin is None) or (fmax is None):
        stf1 = np.fft.fft(np.fft.fftshift(s1[timestamps.min():timestamps.max() + 1]))
        stf2 = np.fft.fft(np.fft.fftshift(s2[timestamps.min():timestamps.max() + 1]))
        nstf = stf1.shape[0]
        sp1 = np.abs(stf1[:int(np.round(nstf / 2.0))]) ** 2
        sp2 = np.abs(stf2[:int(np.round(nstf / 2.0))]) ** 2
        maxsp1 = sp1.max()
        maxsp2 = sp2.max()
        f = np.linspace(0, 0.5, np.round(nstf / 2.0) + 1)[:int(np.round(nstf / 2.0))]
        if fmin is None:
            mask = sp1 > maxsp1 / 100.0
            indmin = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].min()
            fmin = max([0.01, 0.05 * np.floor(f[indmin] / 0.05)])
        if fmax is None:
            mask = sp2 > maxsp2 / 100.0
            indmax = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].max()
            fmax = 0.05 * np.ceil(f[indmax] / 0.05)

    B = fmax - fmin
    R = B / ((fmin + fmax) / 2.0)
    ratio = fmax / fmin
    umax = np.log(ratio)
    teq = nh0 / (fmax * umax)
    if teq > 2 * nh0:
        m0 = (2 * nh0 ** 2) / teq - nh0 + 1
    else:
        m0 = 0
    mu = np.round(nh0 + m0)
    T = 2 * mu - 1

    if n_voices is None:
        nq = np.ceil((B * T * (1 + 2.0 / R) * np.log((1 + R / 2.0) / (1 - R / 2.0))) / 2)
        nmin = nq - nq % 2
        ndflt = 2 ** nextpow2(nmin)
        n_voices = int(ndflt)

    k = np.arange(1, n_voices + 1)
    q = ratio ** (1.0 / (n_voices - 1))
    a = np.exp((k - 1) * np.log(q))
    geo_f = fmin * a

    # Wavelet decomposition computation
    matxte1 = np.zeros((n_voices, tcol), dtype=complex)
    matxte2 = np.zeros((n_voices, tcol), dtype=complex)
    _, _, _, wt1 = scalogram(s1, time_instants=timestamps, waveparams=nh0,
            fmin=fmin, fmax=fmax, n_voices=n_voices)
    _, _, _, wt2 = scalogram(s2, time_instants=timestamps, waveparams=nh0,
            fmin=fmin, fmax=fmax, n_voices=n_voices)
    for ptr in range(n_voices):
        matxte1[ptr, :] = wt1[ptr, :] * np.sqrt(a[n_voices - ptr - 1])
        matxte2[ptr, :] = wt2[ptr, :] * np.sqrt(a[n_voices - ptr - 1])

    umin = -umax
    u = np.linspace(umin, umax, 2 * mu + 1)
    u = u[:(2 * mu)]
    u[mu] = 0
    p = np.arange(2 * n_voices)
    beta = (p / float(n_voices) - 1.0) / (2 * np.log(q))
    l1 = l2 = np.zeros((2 * mu, 2 * n_voices), dtype=complex)
    for m in range(l1.shape[0]):
        l1[m, :] = np.exp(-2 * np.pi * 1j * beta * np.log(lambdak(u[m], K)))
        l2[m, :] = np.exp(-2 * np.pi * 1j * beta * np.log(lambdak(-u[m], K)))

    # Calculate time smoothing window
    if ng0 == 0:
        G = np.ones((2 * mu))
    else:
        a_t = 3
        sigma_t = ng0 * fmax / np.sqrt(2 * a_t * np.log(10))
        a_u = 2 * np.pi ** 2 * sigma_t ** 2 * umax ** 2 / np.log(10)
        G = np.exp(-(a_u * np.log(10) / mu ** 2) * np.arange(-mu, mu) ** 2)

    waf = np.zeros((2 * mu, n_voices))
    tfr = np.zeros((n_voices, tcol))
    S1 = S2 = np.zeros((2 * n_voices,), dtype=complex)
    mx1 = mx2 = np.zeros((2 * n_voices, 2 * mu))

    for ti in range(tcol):
        S1[:n_voices] = matxte1[:, ti]
        mellin1 = np.fft.fftshift(np.fft.ifft(S1))
        mx1 = l1 * mellin1.reshape(1, mellin1.shape[0]).repeat(2 * mu, 0)
        mx1 = np.fft.fft(mx1, axis=0)
        tx1 = mx1[:n_voices, :].T

        S2[:n_voices] = matxte2[:, ti]
        mellin2 = np.fft.fftshift(np.fft.ifft(S2))
        mx2 = l2 * mellin2.reshape(1, mellin2.shape[0]).repeat(2 * mu, 0)
        mx2 = np.fft.fft(mx2, axis=0)
        tx2 = mx2[:n_voices, :].T
        waf = np.real(tx1 * np.conj(tx2)) * G.reshape(G.shape[0], 1).repeat(n_voices, axis=1)
        tfr[:, ti] = np.sum(waf) * geo_f

    t = timestamps
    f = geo_f

    # Normalization
    sp1 = np.fft.fft(hilbert(s1))
    sp2 = np.fft.fft(hilbert(s2))
    indmin = 1 + np.round(fmin * (xrow - 2))
    indmax = 1 + np.round(fmax * (xrow - 2))
    sp1_ana = sp1[indmin:(indmax + 1)]
    sp2_ana = sp2[indmin:(indmax + 1)]
    xx = np.dot(np.real(sp1_ana), np.real(sp2_ana))
    xx += np.dot(np.imag(sp1_ana), np.imag(sp2_ana))
    tfr = tfr * xx / integrate_2d(tfr, t, f) / n_voices
    return tfr, t, f


def umaxdfla_solve(ratio):
    coeffs = [(1.0 - ratio) / 16, (1.0 + ratio) / 2, 1 - ratio]
    roots = np.roots(coeffs)
    return np.min(np.abs(roots - 0))


def unterberger(signal, timestamps=None, form='A', fmin=None, fmax=None,
                n_voices=None):
    """unterberger

    :param signal:
    :param timestamps:
    :param form:
    :param fmin:
    :param fmax:
    :param n_voices:
    :type signal:
    :type timestamps:
    :type form:
    :type fmin:
    :type fmax:
    :type n_voices:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(xrow)

    tcol = timestamps.shape[0]
    x1 = signal.copy()
    x2 = signal.copy()

    s1 = np.real(x1)
    s2 = np.real(x2)
    m = (xrow + (xrow % 2)) / 2
    t = np.arange(xrow) - m - 1
    mt = xrow
    T = xrow

    if (fmin is None) or (fmax is None):
        stf1 = np.fft.fft(np.fft.fftshift(s1[timestamps.min():timestamps.max() + 1]))
        stf2 = np.fft.fft(np.fft.fftshift(s2[timestamps.min():timestamps.max() + 1]))
        nstf = stf1.shape[0]
        sp1 = np.abs(stf1[:int(np.round(nstf / 2.0))]) ** 2
        sp2 = np.abs(stf2[:int(np.round(nstf / 2.0))]) ** 2
        maxsp1 = sp1.max()
        maxsp2 = sp2.max()
        f = np.linspace(0, 0.5, np.round(nstf / 2.0) + 1)[:int(np.round(nstf / 2.0))]
        if fmin is None:
            mask = sp1 > maxsp1 / 100.0
            indmin = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].min()
            fmin = max([0.01, 0.05 * np.floor(f[indmin] / 0.05)])
        if fmax is None:
            mask = sp2 > maxsp2 / 100.0
            indmax = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].max()
            fmax = 0.05 * np.ceil(f[indmax] / 0.05)

    bw = fmax - fmin
    R = bw / (fmin + fmax) * 2.0
    umaxunt = lambda x: (np.sqrt(1 + (x / 2.0) ** 2) + x / 2.0) / (np.sqrt(1 + (x / 2.0) ** 2) - x / 2.0) - fmax / fmin
    umax = newton(umaxunt, 0)
    teq = m / (fmax * umax)
    if teq < 2 * m:
        m0 = np.round((2 * m ** 2) / teq - m) + 1
        m1 = m + m0
        T = 2 * (m + m0) - 1
    else:
        m0 = 0
        m1 = m
    m1 = int(np.round(m1))

    if n_voices is None:
        nq = np.ceil((bw * T * (1 + 2.0 / R) * np.log((1 + R / 2.0) / (1 - R / 2.0))) / 2)
        nmin = nq - nq % 2
        ndflt = 2 ** nextpow2(nmin)
        n_voices = int(ndflt)

    # Geometric sampling for the analyzed spectrum
    k = np.arange(1, n_voices + 1)
    q = (fmax / fmin) ** (1 / (n_voices - 1.0))
    t = np.arange(1, mt + 1) - m - 1
    geo_f = fmin * np.exp((k - 1) * np.log(q))
    tfmatx = np.exp(-2 * 1j * np.dot(t.reshape(t.shape[0], 1),
                                     geo_f.reshape(1, geo_f.shape[0])) * np.pi)
    S1 = np.dot(s1.reshape(1, s1.shape[0]), tfmatx)
    S2 = np.dot(s2.reshape(1, s2.shape[0]), tfmatx)
    S1 = np.append(S1, np.zeros((n_voices,)))
    S2 = np.append(S2, np.zeros((n_voices,)))

    # Mellin tranform of signal
    p = np.arange(2 * n_voices)
    mellin1 = np.fft.fftshift(np.fft.ifft(S1))
    mellin2 = np.fft.fftshift(np.fft.ifft(S2))
    umin = -umax
    du = np.abs(umax - umin) / (2 * m1)
    u = np.linspace(umin, umax - du, (umax - umin) / du)
    u[m1] = 0
    beta = (p / float(n_voices) - 1) / (2 * np.log(q))

    # Computation of Lambda dilations/compressions
    waf = np.zeros((2 * m1, n_voices), dtype=complex)
    _x = -(2 * 1j * np.pi * beta + 0.5)
    for n in range(1, 2 * m1 + 1):
        _y = np.log(np.sqrt(1 + (u[n - 1] / 2.0) ** 2) - u[n - 1] / 2.0)
        mx1 = np.exp(_x * _y) * mellin1
        _y = np.log(np.sqrt(1 + (u[n - 1] / 2.0) ** 2) + u[n - 1] / 2.0)
        mx2 = np.exp(_x * _y) * mellin2
        fx1 = np.fft.fft(np.fft.fftshift(mx1))[:n_voices]
        fx2 = np.fft.fft(np.fft.fftshift(mx2))[:n_voices]
        waf[n - 1, :] = fx1 * np.conj(fx2)
        if form != "A":
            waf[n - 1, :] *= 1 / np.sqrt(1 + (u[n] / 2.0) ** 2)
    waf = np.vstack((waf[m1:(2 * m1), :], waf[:m1, :]))
    waf *= np.repeat(geo_f.reshape((1, geo_f.shape[0])), 2 * m1, axis=0)
    tffr = np.fft.ifft(waf, axis=0)
    tffr = np.real(np.rot90(np.vstack((tffr[m1:(2 * m1 + 1), :],
                                       tffr[:m1, :])), k=-1))

    # conversion from tff to tf using 1d interpolation
    tfr = np.zeros((n_voices, tcol))
    ts2 = (mt - 1.0) / 2
    gamma = np.linspace(-geo_f[n_voices - 1] * ts2,
                        geo_f[n_voices - 1] * ts2, 2 * m1)
    for i in range(n_voices):
        ind = find(np.logical_and(gamma >= -geo_f[i] * ts2,
                                  gamma <= geo_f[i] * ts2))
        x = gamma[ind]
        y = tffr[i, ind]
        xi = (timestamps - ts2 - 1) * geo_f[i]
        tck = splrep(x, y)
        tfr[i, :] = splev(xi, tck).ravel()
    t = timestamps
    f = geo_f[:n_voices].ravel()

    # Normalization
    SP1 = np.fft.fft(hilbert(s1), axis=0)
    SP2 = np.fft.fft(hilbert(s2), axis=0)
    indmin = 1 + int(np.round(fmin * (tcol - 2)))
    indmax = 1 + int(np.round(fmax * (tcol - 2)))
    sp1_ana = SP1[(indmin - 1):indmax]
    sp2_ana = SP2[(indmin - 1):indmax]
    tfr = tfr * np.dot(sp1_ana.T, sp2_ana) / integrate_2d(tfr, t, f) / n_voices

    return tfr, t, f


def d_flandrin(signal, timestamps=None, fmin=None, fmax=None, n_voices=None):
    """d_flandrin

    :param signal:
    :param timestamps:
    :param fmin:
    :param fmax:
    :param n_voices:
    :type signal:
    :type timestamps:
    :type fmin:
    :type fmax:
    :type n_voices:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(xrow)

    tcol = timestamps.shape[0]
    x1 = signal.copy()
    x2 = signal.copy()

    s1 = np.real(x1)
    s2 = np.real(x2)
    m = (xrow + (xrow % 2)) / 2
    t = np.arange(xrow) - m - 1
    mt = xrow
    T = xrow

    if (fmin is None) or (fmax is None):
        stf1 = np.fft.fft(np.fft.fftshift(s1[timestamps.min():timestamps.max() + 1]))
        stf2 = np.fft.fft(np.fft.fftshift(s2[timestamps.min():timestamps.max() + 1]))
        nstf = stf1.shape[0]
        sp1 = np.abs(stf1[:int(np.round(nstf / 2.0))]) ** 2
        sp2 = np.abs(stf2[:int(np.round(nstf / 2.0))]) ** 2
        maxsp1 = sp1.max()
        maxsp2 = sp2.max()
        f = np.linspace(0, 0.5, np.round(nstf / 2.0) + 1)
        if fmin is None:
            mask = sp1 > maxsp1 / 100.0
            indmin = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].min()
            fmin = max([0.01, 0.05 * np.floor(f[indmin] / 0.05)])
        if fmax is None:
            mask = sp2 > maxsp2 / 100.0
            indmax = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].max()
            fmax = 0.05 * np.ceil(f[indmax] / 0.05)
    bw = fmax - fmin
    R = bw / (fmin + fmax) * 2.0
    umax = umaxdfla_solve(fmax / fmin)
    teq = m / (fmax * umax)
    if teq < 2 * m:
        m0 = np.round((2 * m ** 2) / teq - m) + 1
        T = 2 * (m + m0) - 1
    else:
        m0 = 0
    m1 = m + m0
    m1 = int(np.round(m1))

    if n_voices is None:
        nq = np.ceil((bw * T * (1 + 2.0 / R) * np.log((1 + R / 2.0) / (1 - R / 2.0))) / 2)
        nmin = nq - nq % 2
        ndflt = 2 ** nextpow2(nmin)
        n_voices = int(ndflt)

    # Geometric sampling for the analyzed spectrum
    k = np.arange(1, n_voices + 1)
    q = (fmax / fmin) ** (1 / (n_voices - 1.0))
    t = np.arange(1, mt + 1) - m - 1
    geo_f = fmin * np.exp((k - 1) * np.log(q))
    tfmatx = np.exp(-2 * 1j * np.dot(t.reshape(t.shape[0], 1),
                                     geo_f.reshape(1, geo_f.shape[0])) * np.pi)
    S1 = np.dot(s1.reshape(1, s1.shape[0]), tfmatx)
    S2 = np.dot(s2.reshape(1, s2.shape[0]), tfmatx)
    S1 = np.append(S1, np.zeros((n_voices,)))
    S2 = np.append(S2, np.zeros((n_voices,)))

    # Mellin tranform of signal
    p = np.arange(2 * n_voices)
    mellin1 = np.fft.fftshift(np.fft.ifft(S1))
    mellin2 = np.fft.fftshift(np.fft.ifft(S2))
    umin = -umax
    du = np.abs(umax - umin) / (2 * m1)
    u = np.linspace(umin, umax - du, (umax - umin) / du)
    u[m1] = 0
    beta = (p / float(n_voices) - 1) / (2 * np.log(q))

    # Computation of Lambda dilations/compressions
    waf = np.zeros((2 * m1, n_voices), dtype=complex)
    for n in range(1, 2 * m1 + 1):
        mx1 = np.exp(-(2 * 1j * np.pi * beta + 0.5) * 2 * np.log(1 - u[n - 1] / 4)) * mellin1
        mx2 = np.exp(-(2 * 1j * np.pi * beta + 0.5) * 2 * np.log(1 + u[n - 1] / 4)) * mellin2
        fx1 = np.fft.fft(np.fft.fftshift(mx1))[:n_voices]
        fx2 = np.fft.fft(np.fft.fftshift(mx2))[:n_voices]
        waf[n - 1, :] = fx1 * np.conj(fx2)
    waf = np.vstack((waf[m1:(2 * m1), :], waf[:m1, :]))
    waf *= np.repeat(geo_f.reshape((1, geo_f.shape[0])), 2 * m1, axis=0)
    tffr = np.fft.ifft(waf, axis=0)
    tffr = np.real(np.rot90(np.vstack((tffr[m1:(2 * m1 + 1), :],
                                       tffr[:m1, :])), k=-1))

    # conversion from tff to tf using 1d interpolation
    tfr = np.zeros((n_voices, tcol))
    ts2 = (mt - 1.0) / 2
    gamma = np.linspace(-geo_f[n_voices - 1] * ts2,
                        geo_f[n_voices - 1] * ts2, 2 * m1)
    for i in range(n_voices):
        ind = find(np.logical_and(gamma >= -geo_f[i] * ts2,
                                  gamma <= geo_f[i] * ts2))
        x = gamma[ind]
        y = tffr[i, ind]
        xi = (timestamps - ts2 - 1) * geo_f[i]
        tck = splrep(x, y)
        tfr[i, :] = splev(xi, tck).ravel()
    t = timestamps
    f = geo_f.ravel()

    # Normalization
    SP1 = np.fft.fft(hilbert(s1), axis=0)
    SP2 = np.fft.fft(hilbert(s2), axis=0)
    indmin = 1 + int(np.round(fmin * (tcol - 2)))
    indmax = 1 + int(np.round(fmax * (tcol - 2)))
    sp1_ana = SP1[(indmin - 1):indmax]
    sp2_ana = SP2[(indmin - 1):indmax]
    tfr = tfr * np.dot(sp1_ana.T, sp2_ana) / integrate_2d(tfr, t, f) / n_voices

    return tfr, t, f


def bertrand(signal, timestamps=None, fmin=None, fmax=None, n_voices=None):
    """bertrand

    :param signal:
    :param timestamps:
    :param fmin:
    :param fmax:
    :param n_voices:
    :type signal:
    :type timestamps:
    :type fmin:
    :type fmax:
    :type n_voices:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(xrow)

    tcol = timestamps.shape[0]
    x1 = signal.copy()
    x2 = signal.copy()

    s1 = np.real(x1)
    s2 = np.real(x2)
    m = (xrow + (xrow % 2)) / 2
    t = np.arange(xrow) - m - 1
    tmin = 1
    tmax = xrow
    T = tmax - tmin
    mt = xrow

    if (fmin is None) or (fmax is None):
        stf1 = np.fft.fft(np.fft.fftshift(s1[timestamps.min():timestamps.max() + 1]))
        stf2 = np.fft.fft(np.fft.fftshift(s2[timestamps.min():timestamps.max() + 1]))
        nstf = stf1.shape[0]
        sp1 = np.abs(stf1[:int(np.round(nstf / 2.0))]) ** 2
        sp2 = np.abs(stf2[:int(np.round(nstf / 2.0))]) ** 2
        maxsp1 = sp1.max()
        maxsp2 = sp2.max()
        f = np.linspace(0, 0.5, np.round(nstf / 2.0) + 1)
        if fmin is None:
            mask = sp1 > maxsp1 / 100.0
            indmin = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].min()
            fmin = max([0.01, 0.05 * np.floor(f[indmin] / 0.05)])
        if fmax is None:
            mask = sp2 > maxsp2 / 100.0
            indmax = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].max()
            fmax = 0.05 * np.ceil(f[indmax] / 0.05)
    bw = fmax - fmin
    R = bw / (fmin + fmax) * 2.0
    umaxbert = lambda x: np.exp(x) - fmax / fmin
    umax = brenth(umaxbert, 0, 4)
    teq = m / (fmax * umax)
    if teq < mt:
        m0 = np.round((2 * m ** 2) / teq - m) + 1
        m1 = m + m0
        T = 2 * m1 - 1
    else:
        m0 = 0
        m1 = m

    if n_voices is None:
        nq = np.ceil((bw * T * (1 + 2.0 / R) * np.log((1 + R / 2.0) / (1 - R / 2.0))) / 2)
        nmin = nq - nq % 2
        ndflt = 2 ** nextpow2(nmin)
        n_voices = int(ndflt)

    # Geometric sampling for the analyzed spectrum
    k = np.arange(1, n_voices + 1)
    q = (fmax / fmin) ** (1 / (n_voices - 1.0))
    t = np.arange(1, mt + 1) - m - 1
    geo_f = fmin * np.exp((k - 1) * np.log(q))
    tfmatx = np.exp(-2 * 1j * np.dot(t.reshape(t.shape[0], 1),
                                     geo_f.reshape(1, geo_f.shape[0])) * np.pi)
    S1 = np.dot(s1.reshape(1, s1.shape[0]), tfmatx)
    S2 = np.dot(s2.reshape(1, s2.shape[0]), tfmatx)
    S1 = np.append(S1, np.zeros((n_voices,)))
    S2 = np.append(S2, np.zeros((n_voices,)))

    # Mellin tranform of signal
    p = np.arange(2 * n_voices)
    mellin1 = np.fft.fftshift(np.fft.ifft(S1))
    mellin2 = np.fft.fftshift(np.fft.ifft(S2))
    umin = -umax
    du = np.abs(umax - umin) / (2 * m1)
    u = np.linspace(umin, umax - du, (umax - umin) / du)
    u[m1] = 0
    beta = (p / float(n_voices) - 1) / (2 * np.log(q))

    # Computation of P0(t. f, f)
    waf = np.zeros((2 * m1, n_voices), dtype=complex)
    for n in np.hstack((np.arange(1, m1 + 1), np.arange(m1 + 2, 2 * m1 + 1))):
        mx1 = np.exp((-2 * 1j * np.pi * beta + 0.5) * np.log((u[n - 1] / 2) *
            np.exp(-u[n - 1] / 2.0) / np.sinh(u[n - 1] / 2))) * mellin1
        mx2 = np.exp((-2 * 1j * np.pi * beta + 0.5) * np.log((u[n - 1] / 2) *
            np.exp(u[n - 1] / 2.0) / np.sinh(u[n - 1] / 2))) * mellin2
        fx1 = np.fft.fft(np.fft.fftshift(mx1))[:n_voices]
        fx2 = np.fft.fft(np.fft.fftshift(mx2))[:n_voices]
        waf[n - 1, :] = fx1 * np.conj(fx2)
    waf[m1, :] = S1[:n_voices] * np.conj(S2[:n_voices])
    waf = np.vstack((waf[m1:(2 * m1), :], waf[:m1, :]))
    waf *= np.repeat(geo_f.reshape((1, geo_f.shape[0])), 2 * m1, axis=0)
    tffr = np.fft.ifft(waf, axis=0)
    tffr = np.real(np.rot90(np.vstack((tffr[m1:(2 * m1 + 1), :],
                                       tffr[:m1, :])), k=-1))
    # conversion from tff to tf using 1d interpolation
    tfr = np.zeros((n_voices, tcol))
    ts2 = (mt - 1.0) / 2
    gamma = np.linspace(-geo_f[n_voices - 1] * ts2,
                        geo_f[n_voices - 1] * ts2, 2 * m1)
    for i in range(n_voices):
        ind = find(np.logical_and(gamma >= -geo_f[i] * ts2,
                                  gamma <= geo_f[i] * ts2))
        x = gamma[ind]
        y = tffr[i, ind]
        xi = (timestamps - ts2 - 1) * geo_f[i]
        tck = splrep(x, y)
        tfr[i, :] = splev(xi, tck).ravel()
    t = timestamps
    f = geo_f.ravel()

    # Normalization
    SP1 = np.fft.fft(hilbert(s1), axis=0)
    SP2 = np.fft.fft(hilbert(s2), axis=0)
    indmin = 1 + int(np.round(fmin * (tcol - 2)))
    indmax = 1 + int(np.round(fmax * (tcol - 2)))
    sp1_ana = SP1[(indmin - 1):indmax]
    sp2_ana = SP2[(indmin - 1):indmax]

    tfr = tfr * np.dot(sp1_ana.T, sp2_ana) / integrate_2d(tfr, t, f) / n_voices
    return tfr, t, f


def scalogram(signal, fmin=None, fmax=None, n_voices=None, time_instants=None,
              waveparams=None):
    """scalogram

    :param signal:
    :param fmin:
    :param fmax:
    :param n_voices:
    :param time_instants:
    :param waveparams:
    :type signal:
    :type fmin:
    :type fmax:
    :type n_voices:
    :type time_instants:
    :type waveparams:
:return:
:rtype:
    """
    # FIXME: Output from the MATLAB implementation differs significantly.
    if time_instants is None:
        time_instants = np.arange(signal.shape[0])
    if waveparams is None:
        waveparams = np.sqrt(signal.shape[0])
    if n_voices is None:
        n_voices = signal.shape[0]

    s_centered = np.real(signal) - np.real(signal).mean()
    z = hilbert(s_centered)

    if (fmin is None) or (fmax is None):
        stf = np.fft.fft(np.fft.fftshift(z[time_instants.min():time_instants.max() + 1]))
        nstf = stf.shape[0]
        sp = np.abs(stf[:int(np.round(nstf / 2.0))]) ** 2
        maxsp = sp.max()
        f = np.linspace(0, 0.5, np.round(nstf / 2.0) + 1)
        if fmin is None:
            mask = sp > maxsp / 100.0
            indmin = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].min()
            fmin = max([0.01, 0.05 * np.floor(f[indmin] / 0.05)])
        if fmax is None:
            mask = sp > maxsp / 100.0
            indmax = np.arange(mask.shape[0], dtype=int)[mask.astype(bool)].max()
            fmax = 0.05 * np.ceil(f[indmax] / 0.05)

    f = np.logspace(np.log10(fmin), np.log10(fmax), n_voices)
    a = np.logspace(np.log10(fmax / float(fmin)), np.log10(1), n_voices)
    wt = np.zeros((n_voices, time_instants.shape[0]), dtype=complex)
    tfr = np.zeros((n_voices, time_instants.shape[0]), dtype=complex)

    if waveparams > 0:
        for ptr in range(n_voices):
            nha = waveparams * a[ptr]
            tha = np.arange(-np.round(nha), np.round(nha) + 1)
            x = np.exp(-(2 * np.log(10) / (nha ** 2)) * tha ** 2)
            y = np.exp(1j * 2 * np.pi * f[ptr] * tha)
            ha = x * y
            detail = np.convolve(z, ha) / np.sqrt(a[ptr])
            ix = np.arange(round(nha), detail.shape[0] - np.round(nha) + 1,
                           dtype=int)
            wt[ptr, :] = detail[time_instants]
            detail = detail[ix]
            tfr[ptr, :] = detail[time_instants] * np.conj(detail[time_instants])
    elif waveparams == 0:
        for ptr in range(n_voices):
            ha = mexhat(f[ptr])
            nha = (ha.shape[0] - 1) / 2
            detail = np.convolve(z, ha) / np.sqrt(a[ptr])
            ix = np.arange(round(nha) + 1, detail.shape[0] - np.round(nha) + 1)
            detail = detail[ix]
            wt[ptr, :] = detail[time_instants]
            tfr[ptr, :] = detail[time_instants] * np.conj(detail[time_instants])
    elif isinstance(waveparams, np.ndarray):
        rwav, cwav = waveparams.shape
        if cwav > rwav:
            waveparams = waveparams.T
        wavef = np.fft.fft(waveparams, axis=0)
        nwave = waveparams.shape[0]
        f0 = wavef[np.abs(wavef[:nwave / 2]) == np.amax(np.abs(wavef[:nwave / 2]))]
        f0 = ((f0 - 1) * (1 / nwave)).mean()
        a = np.logspace(np.log10(f0 / float(fmin)), np.log10(f0 / float(fmax)), n_voices)
        B = 0.99
        R = B / (1.001 / 2)
        nscale = np.max([128, np.round((B * nwave * (1 + 2.0 / R) * np.log((1 +
            R / 2.0) / (1 - R / 2.0))) / 2)])
        wts = scale(waveparams, a, fmin, fmax, nscale)
        for ptr in range(n_voices):
            ha = wts[ptr, :]
            nha = ha.shape[0] / 2
            detail = np.convolve(z, ha) / np.sqrt(a[ptr])
            detail = detail[int(np.floor(nha)):(detail.shape[0] - np.round(nha))]
            wt[ptr, :] = detail[time_instants]
            tfr[ptr, :] = detail[time_instants] * np.conj(detail[time_instants])

    t = time_instants
    f = f.T
    # Normalization
    SP = np.fft.fft(z, axis=0)
    indmin = 1 + np.round(fmin * (signal.shape[0] - 2))
    indmax = 1 + np.round(fmax * (signal.shape[0] - 2))
    SPana = SP[indmin:(indmax + 1)]
    tfr = np.real(tfr)
    tfr = tfr * (np.linalg.norm(SPana) ** 2) / integrate_2d(tfr, t, f) / n_voices
    return tfr, t, f, wt

if __name__ == '__main__':
    from tftb.generators import altes
    import matplotlib.pyplot as plt
    sig = altes(64, 0.1, 0.45)
    tfr, t, f = smoothed_pseudo_wigner(sig)
    tfr = np.abs(tfr) ** 2
    threshold = np.amax(tfr) * 0.05
    tfr[tfr <= threshold] = 0.0
    plt.imshow(tfr, aspect='auto', origin='bottomleft', extent=[0, 64, 0, 0.5])
    plt.show()
