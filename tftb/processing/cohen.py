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
from tftb.processing.linear import ShortTimeFourierTransform
from tftb.processing.base import BaseTFRepresentation


class Spectrogram(ShortTimeFourierTransform):

    name = "spectrogram"

    def run(self):
        lh = (self.fwindow.shape[0] - 1) // 2
        rangemin = min([round(self.n_fbins / 2.0) - 1, lh])
        starts = -np.min(np.c_[rangemin * np.ones(self.ts.shape), self.ts - 1],
                         axis=1).astype(int)
        ends = np.min(np.c_[rangemin * np.ones(self.ts.shape),
                      self.signal.shape[0] - self.ts], axis=1).astype(int)
        conj_fwindow = np.conj(self.fwindow)
        for icol in range(self.tfr.shape[1]):
            ti = self.ts[icol]
            start = starts[icol]
            end = ends[icol]
            tau = np.arange(start, end + 1).astype(int)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins)
            self.tfr[indices.astype(int), icol] = self.signal[ti + tau - 1] * \
                conj_fwindow[lh + tau] / np.linalg.norm(self.fwindow[lh + tau])
        self.tfr = np.abs(np.fft.fft(self.tfr, axis=0)) ** 2
        return self.tfr, self.ts, self.freqs

    def plot(self, kind='cmap', **kwargs):
        thresh = kwargs.pop("threshold", 0.0)
        super(Spectrogram, self).plot(kind=kind, sqmod=False, threshold=thresh,
                                      **kwargs)


class PageRepresentation(BaseTFRepresentation):

    name = "page representation"

    def run(self):
        for icol in range(self.ts.shape[0]):
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

    name = "pseudo page"

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
        for icol in range(self.ts.shape[0]):
            tau = np.arange(min([self.n_fbins - 1, lh, icol - 1]) + 1)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins) + 1
            self.tfr[indices, icol] = self.fwindow[lh + tau] * \
                self.signal[icol] * np.conj(self.signal[icol - tau])
        self.tfr = np.real(np.fft.fft(self.tfr, axis=0))
        return self.tfr, self.ts, self.freqs


class MargenauHillDistribution(BaseTFRepresentation):

    name = "margenau-hill"

    def run(self):
        for icol in range(self.ts.shape[0]):
            ti = self.ts[icol]
            tau = np.arange(-min((self.n_fbins - ti,
                                  self.signal.shape[0] - ti)) + 1, ti)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins)
            self.tfr[indices, icol] = self.signal[ti] * \
                np.conj(self.signal[ti - tau])

        self.tfr = np.real(np.fft.fft(self.tfr, axis=0))
        return self.tfr, self.ts, self.freqs

    def plot(self, kind='cmap', threshold=0.05, sqmod=True, **kwargs):
        self.tfr = self.tfr[:(self.tfr.shape[0] // 2), :]
        if sqmod:
            self.tfr = np.abs(self.tfr) ** 2
        _threshold = np.amax(self.tfr) * threshold
        self.tfr[self.tfr <= _threshold] = 0.0
        extent = [0, self.ts.max(), 0, 0.5]
        super(MargenauHillDistribution, self).plot(kind=kind, extent=extent, **kwargs)


class PseudoMargenauHillDistribution(MargenauHillDistribution):

    name = "pseudo margenau-hill"

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
        xrow = self.signal.shape[0]
        for icol in range(self.ts.shape[0]):
            start = min([np.round(self.n_fbins / 2.0) - 1, lh, xrow - icol])
            end = min([np.round(self.n_fbins / 2.0) - 1, lh, icol - 1])
            tau = np.arange(-start, end + 1)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins)
            self.tfr[indices, icol] = self.fwindow[lh + tau] * self.signal[icol] * \
                np.conj(self.signal[icol - tau - 1])
        self.tfr = np.fft.fft(self.tfr, axis=0)
        return self.tfr, self.ts, self.freqs


class WignerVilleDistribution(BaseTFRepresentation):

    name = "wigner-ville"

    def run(self):
        tausec = round(self.n_fbins / 2.0)
        winlength = tausec - 1
        taulens = np.min(np.c_[np.arange(self.signal.shape[0]),
                               self.signal.shape[0] - np.arange(self.signal.shape[0]) - 1,
                         winlength * np.ones(self.ts.shape)], axis=1)
        conj_signal = np.conj(self.signal)
        for icol in range(self.ts.shape[0]):
            taumax = taulens[icol]
            tau = np.arange(-taumax, taumax + 1).astype(int)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins).astype(int)
            self.tfr[indices, icol] = self.signal[icol + tau] * \
                conj_signal[icol - tau]
            if (icol <= self.signal.shape[0] - tausec) and (icol >= tausec + 1):
                self.tfr[tausec, icol] = (self.signal[icol + tausec, 0] *
                                          np.conj(self.signal[icol - tausec, 0])) + \
                    (self.signal[icol - tausec, 0] * conj_signal[icol + tausec, 0])
        self.tfr = np.fft.fft(self.tfr, axis=0)
        self.tfr = np.real(self.tfr)
        self.freqs = 0.5 * np.arange(self.n_fbins, dtype=float) / self.n_fbins
        return self.tfr, self.ts, self.freqs

    def plot(self, kind='cmap', threshold=0.05, sqmod=False, **kwargs):
        scale = kwargs.pop("scale", "linear")
        if scale == "log":
            maxi = np.amax(self.tfr)
            mini = max(np.amin(self.tfr), maxi * threshold)
            levels = np.logspace(np.log10(mini), np.log10(maxi), 65)
            kwargs['levels'] = levels
        if sqmod:
            self.tfr = np.abs(self.tfr) ** 2
        else:
            self.tfr = np.abs(self.tfr)
        _threshold = np.amax(self.tfr) * threshold
        self.tfr[self.tfr <= _threshold] = 0.0
        super(WignerVilleDistribution, self).plot(kind=kind, threshold=threshold,
                                                  **kwargs)


class PseudoWignerVilleDistribution(WignerVilleDistribution):

    name = "pseudo winger-ville"

    def run(self):
        lh = (self.fwindow.shape[0] - 1) // 2
        for icol in range(self.ts.shape[0]):
            taumaxvals = (icol, self.signal.shape[0] - icol - 1,
                          np.round(self.n_fbins / 2.0), lh)
            taumax = np.min(taumaxvals)
            tau = np.arange(-taumax, taumax + 1).astype(int)
            indices = np.remainder(self.n_fbins + tau, self.n_fbins).astype(int)
            self.tfr[indices, icol] = self.fwindow[lh + tau] * self.signal[icol + tau] * \
                np.conj(self.signal[icol - tau])
            tau = np.round(self.n_fbins / 2.0)
            if (icol <= self.signal.shape[0] - tau) and (icol >= tau + 1) and (tau <= lh):
                self.tfr[int(tau), icol] = self.fwindow[lh + tau] * \
                    self.signal[icol + tau, 0] * np.conj(self.signal[icol - tau, 0]) + \
                    self.fwindow[lh - tau] * self.signal[icol - tau, 0] * \
                    np.conj(self.signal[icol + tau, 0])
                self.tfr[int(tau), icol] *= 0.5

        self.tfr = np.fft.fft(self.tfr, axis=0)
        return np.real(self.tfr), self.ts, self.freqs

    def plot(self, **kwargs):
        if "kind" not in kwargs:
            kwargs['kind'] = "contour"
        super(PseudoWignerVilleDistribution, self).plot(**kwargs)


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
    if timestamps is None:
        timestamps = np.arange(signal.shape[0])
    if freq_bins is None:
        freq_bins = signal.shape[0]

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
    lg = (twindow.shape[0] - 1) // 2
    lh = (fwindow.shape[0] - 1) // 2
    for icol in range(timestamps.shape[0]):
        ti = timestamps[icol]
        taumax = min([ti + lg - 1, signal.shape[0] - ti + lg,
                      np.round(freq_bins / 2.0) - 1, lh])
        points = np.arange(-min([lg, signal.shape[0] - ti]),
                           min([lg, ti - 1]) + 1)
        g2 = twindow[lg + points]
        g2 = g2 / np.sum(g2)
        tfr[0, icol] = np.sum(g2 * signal[ti - points - 1] * np.conj(signal[ti - points - 1]))
        for tau in range(int(taumax)):
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


if __name__ == '__main__':
    from tftb.generators import anapulse
    sig = anapulse(128)
    t = np.linspace(0, 1, 128)
    spec = WignerVilleDistribution(sig, timestamps=t)
    spec.run()
    spec.plot(kind="contour", scale="log")
