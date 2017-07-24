#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Reassigned TF processing."""

import numpy as np
import scipy.signal as ssig
from tftb.processing.utils import derive_window


def pseudo_wigner_ville(signal, timestamps=None, n_fbins=None, fwindow=None):
    """pseudo_wigner_ville

    :param signal:
    :param timestamps:
    :param n_fbins:
    :param fwindow:
    :type signal:
    :type timestamps:
    :type n_fbins:
    :type fwindow:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(signal.shape[0])
    if n_fbins is None:
        n_fbins = signal.shape[0]
    tcol = timestamps.shape[0]

    if fwindow is None:
        hlength = np.floor(n_fbins / 4.0)
        if hlength % 2 == 0:
            hlength += 1
        fwindow = ssig.hamming(int(hlength))
    elif fwindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing fwindow must have an odd length.')
    lh = (fwindow.shape[0] - 1) // 2
    fwindow = fwindow / fwindow[lh]

    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    tf2 = np.zeros((n_fbins, tcol), dtype=complex)
    dh = derive_window(fwindow)
    for icol in range(tcol):
        ti = timestamps[icol]
        taumax = min([ti - 1, xrow - ti, np.round(n_fbins / 2.0) - 1, lh])
        tau = np.arange(-taumax, taumax + 1)
        indices = np.remainder(n_fbins + tau, n_fbins) + 1
        tfr[indices - 1, icol] = (fwindow[lh + tau] * signal[ti + tau - 1] *
                                  np.conj(signal[ti - tau - 1]))
        tf2[indices - 1, icol] = (dh[lh + tau] * signal[ti + tau - 1] *
                                  np.conj(signal[ti - tau - 1]))
        tau = np.round(n_fbins / 2)
        if (ti <= (xrow - tau)) and (ti > (tau + 1)) and (tau <= lh):
            _x = fwindow[lh + 1 + tau] * signal[ti + tau] * np.conj(signal[ti - tau])
            _y = fwindow[lh + 1 - tau] * signal[ti - tau] * np.conj(signal[ti + tau])
            tfr[tau + 1, icol] = (_x + _y) * 0.5
            _x = dh[lh + 1 + tau] * signal[ti + tau] * np.conj(signal[ti - tau])
            _y = dh[lh + 1 - tau] * signal[ti - tau] * np.conj(signal[ti + tau])
            tf2[tau + 1, icol] = (_x + _y) * 0.5
    tfr = np.real(np.fft.fft(tfr, axis=0))
    tf2 = np.imag(np.fft.fft(tf2, axis=0))
    tfr = tfr.ravel()
    tf2 = tf2.ravel()
    no_warn_mask = tfr != 0
    tf2[no_warn_mask] *= n_fbins / tfr[no_warn_mask] / (2 * np.pi)
    tf2[no_warn_mask] = np.round(tf2[no_warn_mask])
    tfr = tfr.reshape(n_fbins, tcol)
    tf2 = tf2.reshape(n_fbins, tcol)

    rtfr = np.zeros((n_fbins, tcol), dtype=complex)
    tmin = timestamps.min()
    tmax = timestamps.max()
    threshold = 1.0e-6 * (np.abs(signal[tmin:(tmax + 1)]) ** 2).mean()

    for icol in range(tcol):
        for jcol in range(n_fbins):
            if np.abs(tfr[jcol, icol]) > threshold:
                jcolhat = jcol - tf2[jcol, icol]
                jcolhat = int(np.remainder(np.remainder(jcolhat - 1, n_fbins) +
                                           n_fbins, n_fbins))
                jcolhat += 1
                rtfr[jcolhat - 1, icol] += tfr[jcol, icol]
                tf2[jcol, icol] = jcolhat
            else:
                tf2[jcol, icol] = np.inf
                rtfr[jcol, icol] += tfr[jcol, icol]

    return tfr, rtfr, tf2


def pseudo_margenau_hill(signal, timestamps=None, n_fbins=None, fwindow=None):
    """pseudo_margenau_hill

    :param signal:
    :param timestamps:
    :param n_fbins:
    :param fwindow:
    :type signal:
    :type timestamps:
    :type n_fbins:
    :type fwindow:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(signal.shape[0])
    if n_fbins is None:
        n_fbins = signal.shape[0]
    tcol = timestamps.shape[0]

    if fwindow is None:
        hlength = np.floor(n_fbins / 4.0)
        if hlength % 2 == 0:
            hlength += 1
        fwindow = ssig.hamming(hlength)
    elif fwindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing fwindow must have an odd length.')
    lh = (fwindow.shape[0] - 1) / 2
    fwindow = fwindow / fwindow[lh]

    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    tf2 = np.zeros((n_fbins, tcol), dtype=complex)
    dh = derive_window(fwindow)
    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    for icol in range(tcol):
        ti = timestamps[icol]
        start = min([np.round(n_fbins / 2.0) - 1, lh, xrow - ti])
        end = min([np.round(n_fbins / 2.0) - 1, lh, ti - 1])
        tau = np.arange(-start, end + 1)
        indices = np.remainder(n_fbins + tau, n_fbins)
        tfr[indices, icol] = (fwindow[lh + tau] * signal[ti - 1] *
                              np.conj(signal[ti - tau - 1]))
        tf2[indices, icol] = (dh[lh + tau] * signal[ti - 1] *
                              np.conj(signal[ti - tau - 1]))

    tfr = np.fft.fft(tfr, axis=0)
    tf2 = np.fft.fft(tf2, axis=0)
    tfr = tfr.ravel()
    tf2 = tf2.ravel()
    no_warn_mask = tfr != 0
    tf2[no_warn_mask] *= n_fbins / tfr[no_warn_mask] / (2 * np.pi)
    tf2[no_warn_mask] = np.round(tf2[no_warn_mask])
    tfr = np.real(tfr)
    tf2 = np.imag(tf2)
    tfr = tfr.reshape(n_fbins, tcol)
    tf2 = tf2.reshape(n_fbins, tcol)

    rtfr = np.zeros((n_fbins, tcol), dtype=complex)
    threshold = 1.0e-6 * (np.abs(signal) ** 2).mean()

    for icol in range(tcol):
        for jcol in range(n_fbins):
            if np.abs(tfr[jcol, icol]) > threshold:
                jcolhat = jcol - tf2[jcol, icol]
                jcolhat = np.remainder(np.remainder(jcolhat - 1, n_fbins) +
                                       n_fbins, n_fbins)
                jcolhat += 1
                rtfr[jcolhat - 1, icol] += tfr[jcol, icol]
                tf2[jcol, icol] = jcolhat
            else:
                tf2[jcol, icol] = np.inf
                rtfr[jcol, icol] += tfr[jcol, icol]

    return tfr, rtfr, tf2


def pseudo_page(signal, timestamps=None, n_fbins=None, fwindow=None):
    """pseudo_page

    :param signal:
    :param timestamps:
    :param n_fbins:
    :param fwindow:
    :type signal:
    :type timestamps:
    :type n_fbins:
    :type fwindow:
:return:
:rtype:
    """
    if timestamps is None:
        timestamps = np.arange(signal.shape[0])
    if n_fbins is None:
        n_fbins = signal.shape[0]
    tcol = timestamps.shape[0]

    if fwindow is None:
        hlength = np.floor(n_fbins / 4.0)
        if hlength % 2 == 0:
            hlength += 1
        fwindow = ssig.hamming(hlength)
    elif fwindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing fwindow must have an odd length.')
    lh = (fwindow.shape[0] - 1) / 2
    fwindow = fwindow / fwindow[lh]

    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    tf2 = np.zeros((n_fbins, tcol), dtype=complex)
    dh = derive_window(fwindow)

    for icol in range(tcol):
        tau = np.arange(min([n_fbins - 1, lh, icol - 1]) + 1)
        indices = np.remainder(n_fbins + tau, n_fbins) + 1
        tfr[indices, icol] = (fwindow[lh + tau] * signal[icol] *
                              np.conj(signal[icol - tau]))
        tf2[indices, icol] = (dh[lh + tau] * signal[icol] *
                              np.conj(signal[icol - tau]))
        tf2[0, icol] += signal[icol] * np.conj(signal[icol])

    tfr = np.fft.fft(tfr, axis=0)
    tf2 = np.fft.fft(tf2, axis=0)
    tfr = tfr.ravel()
    tf2 = tf2.ravel()
    no_warn_mask = tfr != 0
    tf2[no_warn_mask] *= n_fbins / tfr[no_warn_mask] / (2 * np.pi)
    tf2[no_warn_mask] = np.round(tf2[no_warn_mask])
    tfr = np.real(tfr)
    tf2 = np.imag(tf2)
    tfr = tfr.reshape(n_fbins, tcol)
    tf2 = tf2.reshape(n_fbins, tcol)

    rtfr = np.zeros((n_fbins, tcol), dtype=complex)
    threshold = 1.0e-6 * (np.abs(signal) ** 2).mean()

    for icol in range(tcol):
        for jcol in range(n_fbins):
            if np.abs(tfr[jcol, icol]) > threshold:
                jcolhat = jcol - tf2[jcol, icol]
                jcolhat = np.remainder(np.remainder(jcolhat - 1, n_fbins) +
                                       n_fbins, n_fbins)
                jcolhat += 1
                rtfr[jcolhat - 1, icol] += tfr[jcol, icol]
                tf2[jcol, icol] = jcolhat
            else:
                tf2[jcol, icol] = np.inf
                rtfr[jcol, icol] += tfr[jcol, icol]

    return tfr, rtfr, tf2


def morlet_scalogram(signal, timestamps=None, n_fbins=None, tbp=0.25):
    """morlet_scalogram

    :param signal:
    :param timestamps:
    :param n_fbins:
    :param tbp:
    :type signal:
    :type timestamps:
    :type n_fbins:
    :type tbp:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(signal.shape[0])
    if n_fbins is None:
        n_fbins = signal.shape[0]
    k = 0.001
    tcol = timestamps.shape[0]
    deltat = timestamps[1:] - timestamps[:-1]
    if deltat.min() != deltat.max():
        raise ValueError("Time instants must be regularly sampled.")
    else:
        dt = deltat.min()

    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    tf2 = np.zeros((n_fbins, tcol), dtype=complex)
    M = np.ceil(tbp * n_fbins * np.sqrt(2 * np.log(1 / k)))
    tau = np.arange(M + int(np.round(n_fbins / 2)) + 1)
    hstar = (np.exp(-(tau / (n_fbins * tbp)) ** 2 / 2.0) *
             np.exp(-1j * 2 * np.pi * tau / n_fbins))
    thstar = tau * hstar

    for m in range(1, int(np.round(n_fbins / 2))):
        factor = np.sqrt(m / (tbp * n_fbins))
        for icol in range(tcol):
            ti = timestamps[icol]
            tau_neg = np.arange(1, min([np.ceil(M / m), ti - 1]) + 1).astype(int)
            tau_pos = np.arange(min([np.ceil(M / m), xrow - ti]) + 1).astype(int)
            # positive frequencies
            tfr[m, icol] = np.dot(hstar[m * tau_pos - 1],
                                  signal[ti + tau_pos - 1])
            tf2[m, icol] = np.dot(thstar[m * tau_pos - 1],
                                  signal[ti + tau_pos - 1])
            if tau_neg.shape[0] > 0:
                tfr[m, icol] += np.dot(np.conj(hstar[tau_neg * m]),
                                       signal[ti - tau_neg])
                tf2[m, icol] -= np.dot(np.conj(thstar[tau_neg * m]),
                                       signal[ti - tau_neg])
            # negative frequencies
            tfr[n_fbins - m, icol] = np.dot(np.conj(hstar[tau_pos * m - 1]),
                                            signal[ti + tau_pos - 1])
            tf2[n_fbins - m, icol] = np.dot(np.conj(thstar[tau_pos * m - 1]),
                                            signal[ti + tau_pos - 1])
            if tau_neg.shape[0] > 0:
                tfr[n_fbins - m, icol] += np.dot(hstar[tau_neg * m],
                                                 signal[ti - tau_neg])
                tf2[n_fbins - m, icol] -= np.dot(thstar[tau_neg * m],
                                                 signal[ti - tau_neg])
        tfr[m, :] *= factor
        tf2[m, :] *= factor / m
        tfr[n_fbins - m, :] *= factor
        tf2[n_fbins - m, :] *= factor / m

    m = int(np.round(n_fbins / 2.0))
    factor = np.sqrt(m / (tbp * n_fbins))
    for icol in range(tcol):
        ti = timestamps[icol]
        tau_neg = np.arange(1, min([np.ceil(M / m), ti - 1]) + 1).astype(int)
        tau_pos = np.arange(min([np.ceil(M / m), xrow - ti]) + 1).astype(int)
        tau_pos -= 1
        tau_neg -= 1

        tfr[m, icol] = np.dot(hstar[m * tau_pos], signal[ti + tau_pos])
        tf2[m, icol] = np.dot(thstar[m * tau_pos], signal[ti + tau_pos])
        if tau_neg.shape[0] > 0:
            tfr[m, icol] += np.dot(np.conj(hstar[tau_neg * m]),
                                   signal[ti - tau_neg])
            tf2[m, icol] -= np.dot(np.conj(thstar[tau_neg * m]),
                                   signal[ti - tau_neg])
    tfr[m, :] *= factor
    tf2[m, :] *= factor / m

    tfr, tf2 = tfr.ravel(), tf2.ravel()
    no_warn_mask = tfr != 0
    tf2[no_warn_mask] = tf2[no_warn_mask] / tfr[no_warn_mask]
    tfr = np.abs(tfr) ** 2
    tfr = tfr.reshape(n_fbins, tcol)
    tf2 = tf2.reshape(n_fbins, tcol).astype(complex)

    rtfr = np.zeros((n_fbins, tcol), dtype=complex)
    ex = np.mean(np.abs(signal) ** 2)
    threshold = ex * 1.0e-6
    factor = 2 * np.pi * n_fbins * (tbp ** 2)
    for icol in range(tcol):
        for jcol in range(n_fbins):
            if tfr[jcol, icol] > threshold:
                icolhat = icol + np.round(np.real(tf2[jcol, icol] / dt))
                icolhat = min([max([icolhat, 1]), tcol])
                m = np.remainder(jcol + np.round(n_fbins / 2.0) - 2,
                                 n_fbins)
                m -= np.round(n_fbins / 2.0) + 1
                jcolhat = jcol + np.round(np.imag((m ** 2) * tf2[jcol, icol] / factor))
                jcolhat = (np.remainder(np.remainder(jcolhat - 1, n_fbins) +
                                        n_fbins, n_fbins) + 1)
                rtfr[jcolhat - 1, icolhat - 1] += tfr[jcol, icol]
                tf2[jcol, icol] = jcolhat + 1j * icolhat
            else:
                tf2[jcol, icol] = np.inf * (1 + 1j)
                rtfr[jcol, icol] += tfr[jcol, icol]

    return tfr, rtfr, tf2


def smoothed_pseudo_wigner_ville(signal, timestamps=None, n_fbins=None,
                                 twindow=None, fwindow=None):
    """smoothed_pseudo_wigner_ville

    :param signal:
    :param timestamps:
    :param n_fbins:
    :param twindow:
    :param fwindow:
    :type signal:
    :type timestamps:
    :type n_fbins:
    :type twindow:
    :type fwindow:
:return:
:rtype:
    """
    xrow = signal.shape[0]

    if timestamps is None:
        timestamps = np.arange(signal.shape[0])
    if n_fbins is None:
        n_fbins = signal.shape[0]
    if fwindow is None:
        hlength = np.floor(n_fbins / 4.0)
        hlength += 1 - (hlength % 2)
        fwindow = ssig.hamming(hlength)
    elif fwindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing window must have an odd length.')
    lh = (fwindow.shape[0] - 1) / 2

    if twindow is None:
        glength = np.floor(n_fbins / 4.0)
        glength += 1 - (glength % 2)
        twindow = ssig.hamming(glength)
    elif twindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing window must have an odd length.')
    lg = (twindow.shape[0] - 1) / 2

    tcol = timestamps.shape[0]
    deltat = timestamps[1:] - timestamps[:-1]
    if deltat.min() != deltat.max():
        raise ValueError("Time instants must be regularly sampled.")
    else:
        dt = deltat.min()

    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    tf2 = np.zeros((n_fbins, tcol), dtype=complex)
    tf3 = np.zeros((n_fbins, tcol), dtype=complex)
    dh = derive_window(fwindow)

    for icol in range(tcol):
        ti = timestamps[icol]
        taumax = min([ti + lg - 1, xrow - ti + lg, np.round(n_fbins / 2.0) - 1,
                      lh])
        points = np.arange(-min([lg, xrow - ti]), min([lg, ti - 1]) + 1)
        g2 = twindow[lg + points]
        g2 = g2 / g2.sum()
        tg2 = g2 * points
        xx = signal[ti - 1 - points] * np.conj(signal[ti - 1 - points])
        tfr[0, icol] = (g2 * xx).sum()
        tf2[0, icol] = (tg2 * xx).sum()
        tf3[0, icol] = dh[lh + 1] * tfr[0, icol]

        for tau in range(int(taumax)):
            points = np.arange(-min([lg, xrow - ti - tau]),
                               min([lg, ti - tau - 1]) + 1)
            g2 = twindow[lg + points]
            g2 = g2 / g2.sum()
            tg2 = g2 * points
            xx = signal[ti + tau - 1 - points] * np.conj(signal[ti - tau - 1 - points])
            tfr[tau, icol] = (g2 * xx).sum() * fwindow[lh + tau]
            tf2[tau, icol] = fwindow[lh + tau] * (tg2 * xx).sum()
            tf3[tau, icol] = dh[lh + tau] * (g2 * xx).sum()
            tfr[n_fbins - tau - 1, icol] = (g2 * np.conj(xx)).sum() * fwindow[lh - tau]
            tf2[n_fbins - tau - 1, icol] = (tg2 * np.conj(xx)).sum() * fwindow[lh - tau]
            tf3[n_fbins - tau - 1, icol] = dh[lh - tau] * (g2 * np.conj(xx)).sum()

    tfr = np.real(np.fft.fft(tfr, axis=0)).ravel()
    tf2 = np.real(np.fft.fft(tf2, axis=0)).ravel()
    tf3 = np.imag(np.fft.fft(tf3, axis=0)).ravel()

    no_warn_mask = tfr != 0
    tf2[no_warn_mask] = np.round(tf2[no_warn_mask] / tfr[no_warn_mask] / dt)
    tf3[no_warn_mask] = np.round(n_fbins * tf3[no_warn_mask] /
                                 tfr[no_warn_mask] / (2 * np.pi))
    tfr, tf2, tf3 = [x.reshape(n_fbins, tcol).astype(complex) for x in (tfr, tf2, tf3)]
    tf3 = np.real(tf3)

    rtfr = np.zeros((n_fbins, tcol), dtype=complex)
    ex = np.mean(np.abs(signal) ** 2)
    threshold = ex * 1.0e-6

    for icol in range(tcol):
        for jcol in range(n_fbins):
            if np.abs(tfr[jcol, icol]) > threshold:
                icolhat = min(max([icol - tf2[jcol, icol], 1]), tcol)
                jcolhat = jcol - tf3[jcol, icol]
                jcolhat = (((int(jcolhat) - 1) % n_fbins) + n_fbins) % n_fbins + 1
                rtfr[jcol, icol] += tfr[jcol, icol]
                tf2[jcol, icol] = jcolhat + 1j * icolhat
            else:
                tf2[jcol, icol] = np.inf * (1 + 1j)
                rtfr[jcol, icol] += tfr[jcol, icol]

    return tfr, rtfr, tf2


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
    if window is None:
        wlength = int(np.floor(signal.shape[0] / 4.0))
        wlength += 1 - np.remainder(wlength, 2)
        window = ssig.hamming(wlength)
    elif window.shape[0] % 2 == 0:
        raise ValueError('The smoothing window must have an odd length.')

    tfr = np.zeros((n_fbins, time_samples.shape[0]), dtype=complex)
    tf2 = np.zeros((n_fbins, time_samples.shape[0]), dtype=complex)
    tf3 = np.zeros((n_fbins, time_samples.shape[0]), dtype=complex)
    lh = (window.shape[0] - 1) // 2
    th = window * np.arange(-lh, lh + 1)
    dwin = derive_window(window)

    for icol in range(time_samples.shape[0]):
        ti = time_samples[icol]
        tau = np.arange(-np.min([np.round(n_fbins / 2) - 1, lh, ti]),
                        np.min([np.round(n_fbins / 2) - 1, lh,
                                signal.shape[0] - ti]) + 1).astype(int)
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
    tf3[no_warn_mask] = np.round(np.imag(n_fbins * tf3[no_warn_mask] /
                                         tfr[no_warn_mask] / (2 * np.pi)))

    tfr = np.abs(tfr) ** 2
    tfr = tfr.reshape(n_fbins, time_samples.shape[0])
    tf2 = tf2.reshape(n_fbins, time_samples.shape[0])
    tf3 = tf3.reshape(n_fbins, time_samples.shape[0])
    tf3 = np.real(tf3)

    rtfr = np.zeros((n_fbins, time_samples.shape[0]), dtype=complex)
    ix = np.arange(time_samples.min(), time_samples.max() + 1) - 1
    threshold = 1e-6 * np.mean(np.abs(signal[ix])**2)
    for icol in range(time_samples.shape[0]):
        for jcol in range(n_fbins):
            if np.abs(tfr[jcol, icol]) > threshold:
                icolhat = icol + tf2[jcol, icol]
                icolhat = np.min([np.max([icolhat, 1]), time_samples.shape[0]])
                jcolhat = jcol - tf3[jcol, icol]
                jcolhat = np.remainder(np.remainder(jcolhat - 1, n_fbins) +
                                       n_fbins, n_fbins)
                rtfr[int(jcolhat), int(icolhat) - 1] += tfr[jcol, icol]
                tf2[jcol, icol] = jcolhat + 1j * icolhat
            else:
                tf2[jcol, icol] = np.inf
                rtfr[jcol, icol] += tfr[jcol, icol]
    return tfr, rtfr, tf2


if __name__ == '__main__':
    from tftb.generators import fmlin
    import matplotlib.pyplot as plt
    ts = np.arange(128, step=2)
    sig = fmlin(128, 0.1, 0.4)[0]
    fwindow = ssig.kaiser(17, beta=3 * np.pi)
    _, rtfr, _ = pseudo_wigner_ville(sig, timestamps=ts, n_fbins=64,
                                     fwindow=fwindow)
    rtfr = np.abs(rtfr) ** 2
    threshold = np.amax(rtfr) * 0.05
    rtfr[rtfr <= threshold] = 0.0
    plt.imshow(rtfr[:64, :], aspect='auto', origin="bottomleft",
               extent=[0, 128, 0, 0.5])
    plt.show()
