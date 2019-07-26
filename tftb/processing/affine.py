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
from scipy.optimize import brenth, newton
from scipy.interpolate import splrep, splev
from tftb.generators import mexhat, scale
from tftb.processing.utils import integrate_2d
from tftb.processing.base import BaseTFRepresentation
from tftb.utils import nextpow2


class AffineDistribution(BaseTFRepresentation):

    isaffine = True

    def __init__(self, signal, fmin=None, fmax=None, **kwargs):
        super(AffineDistribution, self).__init__(signal, **kwargs)
        if signal.ndim == 2:
            self.kind = "cross"
            self.x1, self.x2 = self.signal[:, 0].copy(), self.signal[:, 1].copy()
        else:
            self.kind = "auto"
            self.x1 = self.x2 = self.signal.copy()
        self.s1 = np.real(self.x1)
        self.s2 = np.real(self.x2)
        self.m = (self.signal.shape[0] + (self.signal.shape[0] % 2)) // 2
        if (fmin is None) or (fmax is None):
            stf1 = np.fft.fft(np.fft.fftshift(self.s1[self.ts.min():self.ts.max() + 1]))
            stf2 = np.fft.fft(np.fft.fftshift(self.s2[self.ts.min():self.ts.max() + 1]))
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
        self.fmax, self.fmin = fmax, fmin
        self.bw = fmax - fmin
        self.R = self.bw / (fmin + fmax) * 2.0

    def _get_nvoices(self):
        nq = np.ceil((self.bw * self.T * (1 + 2.0 / self.R) * np.log((1 + self.R / 2.0) / (1 - self.R / 2.0))) / 2)
        nmin = nq - nq % 2
        ndflt = 2 ** nextpow2(nmin)
        self.n_voices = int(ndflt)

    def _normalize(self):
        SP1 = np.fft.fft(hilbert(self.s1), axis=0)
        SP2 = np.fft.fft(hilbert(self.s2), axis=0)
        indmin = 1 + int(np.round(self.fmin * (self.ts.shape[0] - 2)))
        indmax = 1 + int(np.round(self.fmax * (self.ts.shape[0] - 2)))
        sp1_ana = SP1[(indmin - 1):indmax]
        sp2_ana = SP2[(indmin - 1):indmax]
        return sp1_ana, sp2_ana

    def _geometric_sample(self):
        k = np.arange(1, self.n_voices + 1)
        self.q = (self.fmax / self.fmin) ** (1 / (self.n_voices - 1.0))
        t = np.arange(1, self.signal.shape[0] + 1) - self.m - 1
        geo_f = self.fmin * np.exp((k - 1) * np.log(self.q))
        self.geo_f = geo_f
        tfmatx = np.exp(-2 * 1j * np.dot(t.reshape(t.shape[0], 1),
                                        geo_f.reshape(1, geo_f.shape[0])) * np.pi)
        S1 = np.dot(self.s1.reshape(1, self.s1.shape[0]), tfmatx)
        S2 = np.dot(self.s2.reshape(1, self.s2.shape[0]), tfmatx)
        S1 = np.append(S1, np.zeros((self.n_voices,)))
        S2 = np.append(S2, np.zeros((self.n_voices,)))
        return S1, S2

    def _mellin_transform(self, S1, S2):
        p = np.arange(2 * self.n_voices)
        mellin1 = np.fft.fftshift(np.fft.ifft(S1))
        mellin2 = np.fft.fftshift(np.fft.ifft(S2))
        umin = -self.umax
        du = np.abs(self.umax - umin) / (2 * self.m1)
        u = np.linspace(umin, self.umax - du, (self.umax - umin) / du)
        u[int(self.m1)] = 0
        self.u = u
        beta = (p / float(self.n_voices) - 1) / (2 * np.log(self.q))
        return beta, mellin1, mellin2

    def plot(self, kind="contour", show_tf=True, threshold=0.05, **kwargs):
        _thresh = np.amax(self.tfr) * threshold
        self.tfr[self.tfr <= _thresh] = 0.0
        freq_y = kwargs.pop("freq_y", np.linspace(self.fmin, self.fmax,
                                                  self.signal.shape[0] / 2))

        super(AffineDistribution, self).plot(kind=kind, show_tf=show_tf,
                                             freq_y=freq_y, **kwargs)

    def _get_interpolated_tf(self, tffr):
        tfr = np.zeros((self.n_voices, self.ts.shape[0]))
        ts2 = (self.signal.shape[0] - 1.0) / 2
        gamma = np.linspace(-self.geo_f[self.n_voices - 1] * ts2,
                            self.geo_f[self.n_voices - 1] * ts2, 2 * self.m1)
        for i in range(self.n_voices):
            condition = np.logical_and(gamma >= -self.geo_f[i] * ts2,
                                       gamma <= self.geo_f[i] * ts2)
            ind = np.nonzero(np.ravel(condition))
            x = gamma[ind]
            y = tffr[i, ind]
            xi = (self.ts - ts2) * self.geo_f[i]
            tck = splrep(x, y)
            tfr[i, :] = splev(xi, tck).ravel()
        t = self.ts
        f = self.geo_f.ravel()
        self.freqs = f
        sp1_ana, sp2_ana = self._normalize()
        if self.kind == "auto":
            multiplier = np.linalg.norm(sp1_ana) ** 2
        else:
            multiplier = np.dot(sp1_ana, sp2_ana)
        tfr = tfr * multiplier / integrate_2d(tfr, t, f) / self.n_voices
        self.tfr = tfr
        return tfr, t, f

    def run(self):
        raise NotImplementedError


class Scalogram(AffineDistribution):
    """Morlet Scalogram.
    """

    name = "scalogram"
    isaffine = False

    def __init__(self, signal, fmin=None, fmax=None, n_voices=None,
                 waveparams=None, **kwargs):
        super(Scalogram, self).__init__(signal, fmin=fmin, fmax=fmax)
        if waveparams is None:
            waveparams = np.sqrt(signal.shape[0])
        if n_voices is None:
            n_voices = self.signal.shape[0]
        self.n_voices = n_voices
        self.waveparams = waveparams
        s_centered = np.real(self.signal) - np.real(self.signal).mean()
        self.z = hilbert(s_centered)
        self.tfr = self.tfr.astype(complex)

    def run(self):
        f = np.logspace(np.log10(self.fmin), np.log10(self.fmax), self.n_voices)
        a = np.logspace(np.log10(self.fmax / float(self.fmin)), np.log10(1), self.n_voices)
        wt = np.zeros((self.n_voices, self.ts.shape[0]), dtype=complex)
        if self.waveparams > 0:
            for ptr in range(self.n_voices):
                nha = self.waveparams * a[ptr]
                tha = np.arange(-int(np.round(nha)), int(np.round(nha)) + 1)
                x = np.exp(-(2 * np.log(10) / (nha ** 2)) * tha ** 2)
                y = np.exp(1j * 2 * np.pi * f[ptr] * tha)
                ha = x * y
                detail = np.convolve(self.z, ha) / np.sqrt(a[ptr])
                ix = np.arange(int(np.round(nha)), detail.shape[0] - int(np.round(nha)) + 1,
                            dtype=int)
                wt[ptr, :] = detail[self.ts]
                detail = detail[ix]
                self.tfr[ptr, :] = detail[self.ts] * np.conj(detail[self.ts])
        elif self.waveparams == 0:
            for ptr in range(self.n_voices):
                ha = mexhat(f[ptr])
                nha = (ha.shape[0] - 1) / 2
                detail = np.convolve(self.z, ha) / np.sqrt(a[ptr])
                ix = np.arange(int(np.round(nha)) + 1, detail.shape[0] - int(np.round(nha)) + 1)
                detail = detail[ix]
                wt[ptr, :] = detail[self.ts]
                self.tfr[ptr, :] = detail[self.ts] * np.conj(detail[self.ts])
        elif isinstance(self.waveparams, np.ndarray):
            rwav, cwav = self.waveparams.shape
            if cwav > rwav:
                self.waveparams = self.waveparams.T
            wavef = np.fft.fft(self.waveparams, axis=0)
            nwave = self.waveparams.shape[0]
            f0 = wavef[np.abs(wavef[:nwave / 2]) == np.amax(np.abs(wavef[:nwave / 2]))]
            f0 = ((f0 - 1) * (1 / nwave)).mean()
            a = np.logspace(np.log10(f0 / float(self.fmin)), np.log10(f0 / float(self.fmax)), self.n_voices)
            B = 0.99
            R = B / (1.001 / 2)
            nscale = np.max([128, np.round((B * nwave * (1 + 2.0 / R) * np.log((1 +
                R / 2.0) / (1 - R / 2.0))) / 2)])
            wts = scale(self.waveparams, a, self.fmin, self.fmax, nscale)
            for ptr in range(self.n_voices):
                ha = wts[ptr, :]
                nha = ha.shape[0] / 2
                detail = np.convolve(self.z, ha) / np.sqrt(a[ptr])
                detail = detail[int(np.floor(nha)):(detail.shape[0] - np.round(nha))]
                wt[ptr, :] = detail[self.ts]
                self.tfr[ptr, :] = detail[self.ts] * np.conj(detail[self.ts])

        # Normalization
        SP = np.fft.fft(self.z, axis=0)
        indmin = 1 + int(np.round(self.fmin * (self.signal.shape[0] - 2)))
        indmax = 1 + int(np.round(self.fmax * (self.signal.shape[0] - 2)))
        SPana = SP[indmin:(indmax + 1)]
        self.tfr = np.real(self.tfr)
        self.tfr = self.tfr * (np.linalg.norm(SPana) ** 2) / integrate_2d(self.tfr, self.ts, f) / self.n_voices
        return self.tfr, self.ts, f, wt


class UnterbergerDistribution(AffineDistribution):

    name = "unterberger"

    def __init__(self, signal, form="A", fmin=None, fmax=None, n_voices=None,
            **kwargs):
        self.form = form
        super(UnterbergerDistribution, self).__init__(signal, fmin=fmin,
                fmax=fmax, n_voices=n_voices, **kwargs)
        umaxunt = lambda x: (np.sqrt(1 + (x / 2.0) ** 2) + x / 2.0) / (np.sqrt(1 + (x / 2.0) ** 2) - x / 2.0) - self.fmax / self.fmin
        self.umax = newton(umaxunt, 0)
        self.m = (self.signal.shape[0] + (self.signal.shape[0] % 2)) // 2
        teq = self.m / (self.fmax * self.umax)
        if teq < 2 * self.m:
            m0 = int(np.round((2 * self.m ** 2) / teq - self.m)) + 1
            m1 = self.m + m0
            self.T = 2 * (self.m + m0) - 1
        else:
            m0 = 0
            m1 = self.m
        self.m1 = int(np.round(m1))
        if n_voices is None:
            self._get_nvoices()
        else:
            self.n_voices = n_voices

    def run(self):
        S1, S2 = self._geometric_sample()
        beta, mellin1, mellin2 = self._mellin_transform(S1, S2)
        # Computation for lambda dilations/compressions
        waf = np.zeros((2 * self.m1, self.n_voices), dtype=complex)
        _x = -(2 * 1j * np.pi * beta + 0.5)
        for n in range(1, 2 * self.m1 + 1):
            _y = np.log(np.sqrt(1 + (self.u[n - 1] / 2.0) ** 2) - self.u[n - 1] / 2.0)
            mx1 = np.exp(_x * _y) * mellin1
            _y = np.log(np.sqrt(1 + (self.u[n - 1] / 2.0) ** 2) + self.u[n - 1] / 2.0)
            mx2 = np.exp(_x * _y) * mellin2
            fx1 = np.fft.fft(np.fft.fftshift(mx1))[:self.n_voices]
            fx2 = np.fft.fft(np.fft.fftshift(mx2))[:self.n_voices]
            waf[n - 1, :] = fx1 * np.conj(fx2)
            if self.form != "A":
                waf[n - 1, :] *= 1 / np.sqrt(1 + (self.u[n] / 2.0) ** 2)
        waf = np.vstack((waf[self.m1:(2 * self.m1), :], waf[:self.m1, :]))
        waf *= np.repeat(self.geo_f.reshape((1, self.geo_f.shape[0])), 2 * self.m1, axis=0)
        tffr = np.fft.ifft(waf, axis=0)
        tffr = np.real(np.rot90(np.vstack((tffr[self.m1:(2 * self.m1 + 1), :],
                                        tffr[:self.m1, :])), k=-1))
        # conversion from tff to tf using 1d interpolation
        tfr = np.zeros((self.n_voices, self.ts.shape[0]))
        ts2 = (self.signal.shape[0] - 1.0) / 2
        gamma = np.linspace(-self.geo_f[self.n_voices - 1] * ts2,
                            self.geo_f[self.n_voices - 1] * ts2, 2 * self.m1)
        for i in range(self.n_voices):
            condition = np.logical_and(gamma >= -self.geo_f[i] * ts2,
                                       gamma <= self.geo_f[i] * ts2)
            ind = np.nonzero(np.ravel(condition))
            x = gamma[ind]
            y = tffr[i, ind]
            xi = (self.ts - ts2 - 1) * self.geo_f[i]
            tck = splrep(x, y)
            tfr[i, :] = splev(xi, tck).ravel()
        t = self.ts
        f = self.freqs = self.geo_f[:self.n_voices].ravel()
        sp1_ana, sp2_ana = self._normalize()

        if self.kind == "auto":
            multiplier = np.linalg.norm(sp1_ana) ** 2
        else:
            multiplier = np.dot(sp1_ana, sp2_ana)

        tfr = tfr * multiplier / integrate_2d(tfr, t, f) / self.n_voices
        self.tfr = tfr
        return tfr, t, f


class DFlandrinDistribution(AffineDistribution):

    name = "d-flandrin"

    def __init__(self, signal, fmin=None, fmax=None, n_voices=None, **kwargs):
        super(DFlandrinDistribution, self).__init__(signal, fmin=fmin,
                                                   fmax=fmax, n_voices=n_voices,
                                                   **kwargs)
        self.umax = umaxdfla_solve(self.fmax / self.fmin)
        teq = self.m / (self.fmax * self.umax)
        if teq < 2 * self.m:
            m0 = int(round((2 * self.m ** 2) / teq - self.m)) + 1
            self.T = 2 * (self.m + m0) - 1
        else:
            m0 = 0
        self.m1 = int(np.round(self.m + m0))
        if n_voices is None:
            self._get_nvoices()
        else:
            self.n_voices = n_voices

    def run(self):
        S1, S2 = self._geometric_sample()
        beta, mellin1, mellin2 = self._mellin_transform(S1, S2)
        # Computation of Lambda dilations/compressions
        waf = np.zeros((2 * self.m1, self.n_voices), dtype=complex)
        for n in range(1, 2 * self.m1 + 1):
            mx1 = np.exp(-(2 * 1j * np.pi * beta + 0.5) * 2 * np.log(1 - self.u[n - 1] / 4)) * mellin1
            mx2 = np.exp(-(2 * 1j * np.pi * beta + 0.5) * 2 * np.log(1 + self.u[n - 1] / 4)) * mellin2
            fx1 = np.fft.fft(np.fft.fftshift(mx1))[:self.n_voices]
            fx2 = np.fft.fft(np.fft.fftshift(mx2))[:self.n_voices]
            waf[n - 1, :] = fx1 * np.conj(fx2)
        waf = np.vstack((waf[self.m1:(2 * self.m1), :], waf[:self.m1, :]))
        waf *= np.repeat(self.geo_f.reshape((1, self.geo_f.shape[0])), 2 * self.m1, axis=0)
        tffr = np.fft.ifft(waf, axis=0)
        tffr = np.real(np.rot90(np.vstack((tffr[self.m1:(2 * self.m1 + 1), :],
                                        tffr[:self.m1, :])), k=-1))
        # conversion from tff to tf using 1d interpolation
        return self._get_interpolated_tf(tffr)


class BertrandDistribution(AffineDistribution):

    name = "bertrand"

    def __init__(self, signal, fmin=None, fmax=None, n_voices=None, **kwargs):
        super(BertrandDistribution, self).__init__(signal, fmin=fmin,
                                                   fmax=fmax, n_voices=n_voices,
                                                   **kwargs)
        umaxbert = lambda x: np.exp(x) - self.fmax / self.fmin
        self.umax = brenth(umaxbert, 0, 4)
        teq = self.m / (self.fmax * self.umax)
        if teq < self.signal.shape[0]:
            m0 = int(np.round((2 * self.m ** 2) / teq - self.m)) + 1
            m1 = self.m + m0
            self.T = 2 * m1 - 1
        else:
            m0 = 0
            m1 = self.m
        self.m1 = m1
        if n_voices is None:
            self._get_nvoices()
        else:
            self.n_voices = n_voices

    def run(self):
        S1, S2 = self._geometric_sample()
        beta, mellin1, mellin2 = self._mellin_transform(S1, S2)
        # Computation of P0(t. f, f)
        waf = np.zeros((2 * int(self.m1), self.n_voices), dtype=complex)
        for n in np.hstack((np.arange(1, self.m1 + 1), np.arange(self.m1 + 2, 2 * self.m1 + 1))):
            mx1 = np.exp((-2 * 1j * np.pi * beta + 0.5) * np.log((self.u[n - 1] / 2) *
                np.exp(-self.u[n - 1] / 2.0) / np.sinh(self.u[n - 1] / 2))) * mellin1
            mx2 = np.exp((-2 * 1j * np.pi * beta + 0.5) * np.log((self.u[n - 1] / 2) *
                np.exp(self.u[n - 1] / 2.0) / np.sinh(self.u[n - 1] / 2))) * mellin2
            fx1 = np.fft.fft(np.fft.fftshift(mx1))[:self.n_voices]
            fx2 = np.fft.fft(np.fft.fftshift(mx2))[:self.n_voices]
            waf[n - 1, :] = fx1 * np.conj(fx2)
        waf[self.m1, :] = S1[:self.n_voices] * np.conj(S2[:self.n_voices])
        waf = np.vstack((waf[self.m1:(2 * self.m1), :], waf[:self.m1, :]))
        waf *= np.repeat(self.geo_f.reshape((1, self.geo_f.shape[0])), 2 * self.m1, axis=0)
        tffr = np.fft.ifft(waf, axis=0)
        tffr = np.real(np.rot90(np.vstack((tffr[self.m1:(2 * self.m1 + 1), :],
                                        tffr[:self.m1, :])), k=-1))
        # conversion from tff to tf using 1d interpolation
        return self._get_interpolated_tf(tffr)


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
    _, _, _, wt1 = Scalogram(s1, time_instants=timestamps, waveparams=nh0,
            fmin=fmin, fmax=fmax, n_voices=n_voices).run()
    _, _, _, wt2 = Scalogram(s2, time_instants=timestamps, waveparams=nh0,
            fmin=fmin, fmax=fmax, n_voices=n_voices).run()
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
    indmin = 1 + int(np.round(fmin * (xrow - 2)))
    indmax = 1 + int(np.round(fmax * (xrow - 2)))
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


if __name__ == '__main__':
    from tftb.generators import altes
    import matplotlib.pyplot as plt
    sig = altes(64, 0.1, 0.45)
    tfr, timestamps, frequencies = smoothed_pseudo_wigner(sig)
    tfr = np.abs(tfr) ** 2
    threshold = np.amax(tfr) * 0.05
    tfr[tfr <= threshold] = 0.0
    plt.imshow(tfr, aspect='auto', origin='bottomleft', extent=[0, 64, 0, 0.5])
    plt.show()
