#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Reassigned TF processing."""

import warnings
import numpy as np
import scipy.signal as ssig
from tftb.processing.utils import derive_window
from tftb.utils import nextpow2


def spectrogram(signal, time_samples=None, n_fbins=None, window=None):
    """Compute the spectrogram and reassigned spectrogram.

    :param signal: signal to be analzsed
    :param time_samples: time instants (default: np.arange(len(signal)))
    :param n_fbins: number of frequency bins (default: len(signal))
    :param window: frequency smoothing window (default: Hamming with \
        size=len(signal)/4)
    :type signal: array-like
    :type time_samples: array-like
    :type n_fbins: int
    :type window: array-like
    :return: spectrogram, reassigned specstrogram and matrix of reassignment
    vectors
    :rtype: tuple(array-like)
    """

    if time_samples is None:
        time_samples = np.arange(signal.shape[0])
    elif np.unique(np.diff(time_samples)).shape[0] > 1:
        raise ValueError('Time instants must be regularly sampled.')
    if n_fbins is None:
        n_fbins = signal.shape[0]
    elif 2 ** nextpow2(n_fbins) != n_fbins:
        msg = "For faster computations, n_fbins should be a power of 2."
        warnings.warn(msg, UserWarning)
    if window is None:
        wlength = int(np.floor(signal.shape[0] / 4.0))
        wlength += 1 - np.remainder(wlength, 2)
        window = ssig.hamming(wlength)
    elif window.shape[0] % 2 == 0:
        raise ValueError('The smoothing window must have an odd length.')

    tfr = np.zeros((n_fbins, time_samples.shape[0]))
    tf2 = np.zeros((n_fbins, time_samples.shape[0]))
    tf3 = np.zeros((n_fbins, time_samples.shape[0]))
    lh = (window.shape[0] - 1) / 2
    th = window * np.arange(-lh, lh + 1)
    dwin = derive_window(window)

    for icol in xrange(time_samples.shape[0]):
        ti = time_samples[icol]
        tau = np.arange(-np.min([np.round(n_fbins / 2) - 1, lh, ti - 1]),
                        np.min([np.round(n_fbins / 2) - 1, lh, signal.shape[0] - ti]) + 1)
        indices = np.remainder(n_fbins + tau, n_fbins)
        norm_h = np.linalg.norm(window[lh + tau], ord=2)
        tfr[indices, icol] = signal[ti + tau - 1] * np.conj(window[lh + tau]) / norm_h
        tf2[indices, icol] = signal[ti + tau - 1] * np.conj(th[lh + tau]) / norm_h
        tf3[indices, icol] = signal[ti + tau - 1] * np.conj(dwin[lh + tau]) / norm_h

    tfr = np.fft.fft(tfr, axis=0).ravel()
    tf2 = np.fft.fft(tf2, axis=0).ravel()
    tf3 = np.fft.fft(tf3, axis=0).ravel()

    no_warn_mask = tfr != 0
    tf2[no_warn_mask] = np.round(np.real(tf2[no_warn_mask] / tfr[no_warn_mask]))
    tf3[no_warn_mask] = np.round(np.imag(n_fbins * tf3[no_warn_mask] / tfr[no_warn_mask] / (2 * np.pi)))

    tfr = np.abs(tfr) ** 2
    tfr = tfr.reshape(n_fbins, time_samples.shape[0])
    tf2 = tf2.reshape(n_fbins, time_samples.shape[0])
    tf3 = tf3.reshape(n_fbins, time_samples.shape[0])
    tf2 = np.real(tf2)
    tf3 = np.real(tf3)

    rtfr = np.zeros((n_fbins, time_samples.shape[0]))
    ix = np.arange(time_samples.min(), time_samples.max() + 1) - 1
    threshold = 1e-6 * np.mean(np.abs(signal[ix])**2)
    for icol in xrange(time_samples.shape[0]):
        for jcol in xrange(n_fbins):
            if np.abs(tfr[jcol, icol]) > threshold:
                icolhat = icol + tf2[jcol, icol]
                icolhat = np.min([np.max([icolhat, 1]), time_samples.shape[0]])
                jcolhat = jcol - tf3[jcol, icol]
                jcolhat = np.remainder(np.remainder(jcolhat - 1, n_fbins) + n_fbins,
                                       n_fbins)
                rtfr[int(jcolhat), int(icolhat) - 1] += tfr[jcol, icol]
                tf2[jcol, icol] = jcolhat + 1j * icolhat
            else:
                tf2[jcol, icol] = np.inf
                rtfr[jcol, icol] += tfr[jcol, icol]
    return tfr, rtfr, tf2


if __name__ == '__main__':
    from tftb.generators.api import fmlin
    from scipy.signal import kaiser
    signal = fmlin(128, 0.1, 0.4)[0]
    time_samples = np.arange(1, 129, step=2)
    window = kaiser(17, beta=3 * np.pi)
    a, b, c = spectrogram(signal, time_samples, 64, window)
