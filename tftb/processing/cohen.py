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

import numpy as np
from tftb.utils import init_default_args
from tftb.processing.linear import ShortTimeFourierTransform
from tftb.processing.base import BaseTFRepresentation


class Spectrogram(ShortTimeFourierTransform):

    name = "spectrogram"

    def run(self):
        super(Spectrogram, self).run()
        self.tfr = np.abs(self.tfr) ** 2

    def plot(self, kind='cmap', **kwargs):
        super(Spectrogram, self).plot(kind=kind, sqmod=False, threshold=0,
                                      **kwargs)


class PageRepresentation(BaseTFRepresentation):

    name = "page representation"

    def run(self):
        for icol in xrange(self.ts.shape[0]):
            ti = self.ts[icol]
            tau = np.arange(-min([self.n_fbins - ti,
                                  self.signal.shape[0] - ti]), ti)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins)
            self.tfr[indices, icol] = np.dot(self.signal[ti],
                                             np.conj(self.signal[ti - tau - 1]))
        self.tfr = np.real(np.fft.fft(self.tfr, axis=0))
        return self.tfr, self.ts, self.freqs

    def plot(self, kind='cmap', threshold=0.05, sqmod=True, **kwargs):
        self.tfr = self.tfr[:(self.tfr.shape[0] / 2), :]
        self.tfr = np.abs(self.tfr) ** 2
        _threshold = np.amax(self.tfr) * threshold
        self.tfr[self.tfr <= _threshold] = 0.0
        super(PageRepresentation, self).plot(kind=kind, **kwargs)


class PseudoPageRepresentation(PageRepresentation):

    def _make_window(self):
        hlength = np.floor(self.n_fbins / 4.0)
        if hlength % 2 == 0:
            hlength += 1
        from scipy.signal import hamming
        fwindow = hamming(hlength)
        lh = (fwindow.shape[0] - 1) / 2
        return fwindow / fwindow[lh]

    def run(self):
        lh = (self.fwindow.shape[0] - 1) / 2
        for icol in xrange(self.ts.shape[0]):
            tau = np.arange(min([self.n_fbins - 1, lh, icol - 1]) + 1)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins) + 1
            self.tfr[indices, icol] = self.fwindow[lh + tau] * self.signal[icol] * np.conj(
                    self.signal[icol - tau])
        self.tfr = np.real(np.fft.fft(self.tfr, axis=0))


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
    timestamps, n_fbins = init_default_args(signal, timestamps=timestamps,
                                            n_fbins=n_fbins)
    tcol = timestamps.shape[0]

    if fwindow is None:
        hlength = np.floor(n_fbins / 4.0)
        if hlength % 2 == 0:
            hlength += 1
        from scipy.signal import hamming
        fwindow = hamming(hlength)
    elif fwindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing fwindow must have an odd length.')
    lh = (fwindow.shape[0] - 1) / 2
    fwindow = fwindow / fwindow[lh]

    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    for icol in xrange(tcol):
        start = min([np.round(n_fbins / 2.0) - 1, lh, xrow - icol])
        end = min([np.round(n_fbins / 2.0) - 1, lh, icol - 1])
        tau = np.arange(-start, end + 1)
        indices = np.remainder(n_fbins + tau, n_fbins)
        tfr[indices, icol] = fwindow[lh + tau] * signal[icol] * np.conj(signal[icol - tau - 1])

    tfr = np.real(np.fft.fft(tfr, axis=0))

    if n_fbins % 2 == 0:
        freq = np.hstack((np.arange(n_fbins / 2), np.arange(-n_fbins / 2, 0))) / n_fbins
    else:
        freq = np.hstack((np.arange((n_fbins - 1) / 2), np.arange(-(n_fbins - 1) / 2, 0))) / n_fbins

    return tfr, timestamps, freq


def smoothed_pseudo_wigner_ville(signal, timestamps=None, freq_bins=None,
                                twindow=None, fwindow=None):
    """Smoothed Pseudo Wigner-Ville time-frequency distribution.
    :param signal: signal to be analyzed
    :param timestamps: time instants of the signal
    :param freq_bins: number of frequency bins
    :param twindow: time smoothing window
    :param fwindow: frequency smoothing window
    :type signal: array-like
    :type timestamps: array-like
    :type freq_bins: int
    :type twindow: array-like
    :type fwindow: array-like
    :return: Smoothed pseudo Wigner Ville distribution
    :rtype: array-like
    """
    timestamps, freq_bins = init_default_args(signal, timestamps=timestamps,
                                              n_fbins=freq_bins)

    if fwindow is None:
        winlength = np.floor(freq_bins / 4.0)
        winlength = winlength + 1 - np.remainder(winlength, 2)
        from scipy.signal import hamming
        fwindow = hamming(int(winlength))
    elif fwindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing fwindow must have an odd length.')

    if twindow is None:
        timelength = np.floor(freq_bins / 10.0)
        timelength += 1 - np.remainder(timelength, 2)
        from scipy.signal import hamming
        twindow = hamming(int(timelength))
    elif twindow.shape[0] % 2 == 0:
        raise ValueError('The smoothing fwindow must have an odd length.')

    tfr = np.zeros((freq_bins, timestamps.shape[0]), dtype=complex)
    lg = (twindow.shape[0] - 1) / 2
    lh = (fwindow.shape[0] - 1) / 2
    for icol in xrange(timestamps.shape[0]):
        ti = timestamps[icol]
        taumax = min([ti + lg - 1, signal.shape[0] - ti + lg,
                      np.round(freq_bins / 2.0) - 1, lh])
        points = np.arange(-min([lg, signal.shape[0] - ti]),
                           min([lg, ti - 1]) + 1)
        g2 = twindow[lg + points]
        g2 = g2 / np.sum(g2)
        tfr[0, icol] = np.sum(g2 * signal[ti - points - 1] * np.conj(signal[ti - points - 1]))
        for tau in xrange(int(taumax)):
            points = np.arange(-min([lg, signal.shape[0] - ti - tau]),
                               min([lg, ti - 1 - tau]) + 1)
            g2 = twindow[lg + points]
            g2 = g2 / np.sum(g2)
            R = np.sum(g2 * signal[ti + tau - points - 1] * np.conj(signal[ti - tau - points - 1]))
            tfr[1 + tau, icol] = fwindow[lh + tau + 1] * R
            R = np.sum(g2 * signal[ti - tau - points - 1] * np.conj(signal[ti + tau - points - 1]))
            tfr[freq_bins - tau - 1, icol] = fwindow[lh - tau + 1] * R
        tau = np.round(freq_bins / 2.0)
        if (ti <= signal.shape[0] - tau) and (ti >= tau + 1) and (tau <= lh):
            points = np.arange(-min([lg, signal.shape[0] - ti - tau]),
                               min([lg, ti - 1 - tau]) + 1)
            g2 = twindow[lg + 1 + points]
            g2 = g2 / np.sum(g2)
            _x = np.sum(g2 * signal[ti + tau - points] * np.conj(signal[ti - tau - points]))
            _x *= fwindow[lh + tau + 1]
            _y = np.sum(g2 * signal[ti - tau - points] * np.conj(signal[ti + tau - points]))
            _y *= fwindow[lh - tau + 1]
            tfr[tau, icol] = (_x + _y) * 0.5
    tfr = np.fft.fft(tfr, axis=0)
    return np.real(tfr)


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
    time_samples, freq_bins = init_default_args(signal, timestamps=time_samples,
                                                n_fbins=freq_bins)

    if window is None:
        winlength = np.floor(freq_bins / 4.0)
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
    time_samples, freq_bins = init_default_args(signal, timestamps=time_samples,
                                                n_fbins=freq_bins)

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


def margenau_hill(signal, timestamps=None, n_fbins=None):
    """Margenau-Hill time frequency distribution.

    :param signal: Signal to be analyzed.
    :param timestamps: Time instants
    :param n_fbins: number of frequency bins
    :type signal: array-like
    :type timestamps: array-like
    :type n_fbins: int
    :return: Margenau Hill representation, timestamps, frequency vector
    :rtype: tuple
    """
    timestamps, n_fbins = init_default_args(signal, timestamps=timestamps,
                                            n_fbins=n_fbins)
    xrow = signal.shape[0]
    tcol = timestamps.shape[0]

    tfr = np.zeros((n_fbins, tcol), dtype=complex)
    for icol in xrange(tcol):
        ti = timestamps[icol]
        tau = np.arange(-min((n_fbins - ti, xrow - ti)) + 1, ti)
        indices = np.remainder(n_fbins + tau, n_fbins)
        tfr[indices, icol] = signal[ti] * np.conj(signal[ti - tau])

    tfr = np.real(np.fft.fft(tfr, axis=0))
    if n_fbins % 2 == 0:
        freq = np.hstack((np.arange(n_fbins / 2), np.arange(-n_fbins / 2, 0))) / n_fbins
    else:
        freq = np.hstack((np.arange((n_fbins - 1) / 2), np.arange(-(n_fbins - 1) / 2, 0))) / n_fbins

    return tfr, timestamps, freq


if __name__ == '__main__':
    from tftb.generators import fmlin
    sig = fmlin(128, 0.1, 0.4)[0]
    spec = PseudoPageRepresentation(sig, n_fbins=64)
    spec.run()
    spec.plot()
