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
from scipy.optimize import brenth
from scipy.interpolate import splrep, splev
from tftb.generators.api import mexhat, scale
from tftb.processing.utils import integrate_2d
from tftb.utils import nextpow2


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
    for i in xrange(n_voices):
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
    wt = np.zeros((n_voices, time_instants.shape[0]))
    tfr = np.zeros((n_voices, time_instants.shape[0]), dtype=complex)

    if waveparams > 0:
        for ptr in xrange(n_voices):
            nha = waveparams * a[ptr]
            tha = np.arange(-np.round(nha), np.round(nha) + 1)
            x = np.exp(-(2 * np.log(10) / (nha ** 2)) * tha ** 2)
            y = np.exp(1j * 2 * np.pi * f[ptr] * tha)
            ha = x * y
            detail = np.convolve(z, ha) / np.sqrt(a[ptr])
            ix = np.arange(round(nha), detail.shape[0] - np.round(nha) + 1,
                           dtype=int)
            detail = detail[ix]
            tfr[ptr, :] = detail[time_instants] * np.conj(detail[time_instants])
    elif waveparams == 0:
        for ptr in xrange(n_voices):
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
        for ptr in xrange(n_voices):
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
    from tftb.generators.api import altes
    import matplotlib.pyplot as plt
    sig = altes(64, 0.1, 0.45)
    tfr, t, f = bertrand(sig)
    tfr = np.abs(tfr) ** 2
    threshold = np.amax(tfr) * 0.05
    tfr[tfr <= threshold] = 0.0
    plt.imshow(tfr, aspect='auto', origin='bottomleft', extent=[0, 64, 0, 0.5])
    plt.show()
