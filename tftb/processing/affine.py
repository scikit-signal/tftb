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
from scipy.signal import hilbert
from tftb.generators.api import mexhat, scale
from tftb.processing.utils import integrate_2d


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
    tfr, t, f, wt = scalogram(sig)
    tfr = np.abs(tfr) ** 2
    threshold = np.amax(tfr) * 0.05
    tfr[tfr <= threshold] = 0.0
    plt.imshow(tfr, extent=[0, 64, 0, 0.5], origin='bottomleft', aspect='auto')
    plt.show()
