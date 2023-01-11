#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Linear Time Frequency Processing.
"""

import numpy as np
import scipy.signal
from scipy.signal.windows import hamming

from tftb.processing.base import BaseTFRepresentation
from tftb.utils import nearest_odd, divider, modulo, izak


class ShortTimeFourierTransform(BaseTFRepresentation):
    """Short time Fourier transform."""

    name = "stft"

    def __init__(self, signal, timestamps=None, n_fbins=None, fwindow=None):
        """Create a ShortTimeFourierTransform object.

        :param signal: Signal to be analyzed.
        :param timestamps: Time instants of the signal (default:
            ``np.arange(len(signal))``)
        :param n_fbins: Number of frequency bins (default: ``len(signal)``)
        :param fwindow: Frequency smoothing window (default: Hamming window of
            length ``len(signal) / 4``)
        :type signal: array-like
        :type timestamps: array-like
        :type n_fbins: int
        :type fwindow: array-like
        :return: ShortTimeFourierTransform object
        :Example:

        >>> from tftb.generators import fmconst
        >>> sig = np.r_[fmconst(128, 0.2)[0], fmconst(128, 0.4)[0]]
        >>> stft = ShortTimeFourierTransform(sig)
        >>> tfr, t, f = stft.run()
        >>> stft.plot() #doctest: +SKIP

        .. plot:: docstring_plots/processing/stft.py
        """
        super().__init__(signal=signal,
                         n_fbins=n_fbins,
                         timestamps=timestamps,
                         fwindow=fwindow)
        # super(ShortTimeFourierTransform, self).__init__(x=x,
        #                                                 n_fbins=n_fbins,
        #                                                 timestamps=timestamps,
        #                                                 fwindow=fwindow)

    def run(self, **kwargs):
        r"""Compute the STFT according to:

        .. math:: X[m, w] = \sum_{n=-\infty}^{\infty}x[n]w[n - m]e^{-j\omega n}

        Where :math:`w` is a Hamming window."""
        """
        
        scipy.x.stft(
            x,
            fs=1.0,
            window='hann',
            nperseg=256,
            noverlap=None,
            nfft=None,
            detrend=False,
            return_onesided=True,
            boundary='zeros',
            padded=True,
            axis=- 1,
            scaling='spectrum')
        """

        self.freqs, self.ts, self.tfr = scipy.signal.stft(self.signal, **kwargs)
        return self.tfr, self.ts, self.freqs

    def plot(self, ax=None, kind='cmap', sqmod=True, threshold=0.05, **kwargs):
        """Display the spectrogram of an STFT.

        :param ax: axes object to draw the plot on. If None(default), one will
            be created.
        :param kind: Choice of visualization type, either "cmap"(default) or "contour".
        :param sqmod: Whether to take squared modulus of TFR before plotting.
            (Default: True)
        :param threshold: Percentage of the maximum value of the TFR, below
            which all values are set to zero before plotting.
        :param **kwargs: parameters passed to the plotting function.
        :type ax: matplotlib.axes.Axes object
        :type kind: str
        :type sqmod: bool
        :type threshold: float
        :return: None
        :rtype: None
        """
        self.tfr = self.tfr[:int(self.n_fbins / 2.0), :]
        self.freqs = self.freqs[:int(self.n_fbins / 2.0)]
        if sqmod:
            self.tfr = np.abs(self.tfr) ** 2
        _threshold = np.amax(self.tfr) * threshold
        self.tfr[self.tfr <= _threshold] = 0.0
        # super(ShortTimeFourierTransform, self).plot(ax=ax, kind=kind,
        #                                             threshold=threshold,
        #                                             **kwargs)
        super().plot(ax=ax, kind=kind, threshold=threshold, **kwargs)


def gabor(signal, n_coeff=None, q_oversample=None, window=None):
    """Compute the Gabor representation of a signal.

    :param signal: Singal to be analyzed.
    :param n_coeff: number of Gabor coefficients in time.
    :param q_oversample: Degree of oversampling
    :param window: Synthesis window
    :type signal: array-like
    :type n_coeff: integer
    :type q_oversample: int
    :type window: array-like
    :return: Tuple of Gabor coefficients, biorthogonal window associated with the synthesis window.
    :rtype: tuple
    """
    if n_coeff is None:
        n_coeff, _ = divider(signal.shape[0])
    if q_oversample is None:
        q_oversample, _ = divider(n_coeff)
    if window is None:
        window = np.exp(np.log(0.005) * np.linspace(-1, 1, nearest_odd(n_coeff)) ** 2)
        window = window / np.linalg.norm(window)
    m = int(q_oversample * signal.shape[0] / float(n_coeff))
    mb = int(signal.shape[0] / float(n_coeff))
    nb = int(signal.shape[0] / float(m))

    # Zak transform?
    nh = window.shape[0]
    if nh % 2 == 0:
        raise ValueError("The window function should have an odd length.")
    alpha = np.round((2 * signal.shape[0] / float(n_coeff) - 1 - nh) / (2 * q_oversample))
    hn1 = np.zeros((signal.shape[0],))
    start = np.round(((signal.shape[0] - (nh - 1))) / 2) - alpha
    end = np.round((signal.shape[0] + nh - 1) / 2) - alpha
    hn1[np.arange(start - 1, end).astype(int)] = window

    msig = hn1.reshape(int(nb), int(m), order='F')
    dzth = np.fft.fft(msig.T, axis=0) / np.sqrt(m)
    mzh = np.zeros((m, mb))
    x = np.arange(1, m + 1, dtype=float)
    for l in range(q_oversample):  # NOQA: E741
        mod = modulo(x - l * m / q_oversample, m).astype(int)
        mzh += np.abs(dzth[mod - 1, :]) ** 2

    mzh[mzh < np.spacing(1)] = 1

    # Za transform of biorthogonal dual frame window gam
    dztgam = dzth / mzh
    gam = np.real(izak(dztgam)) / signal.shape[0]

    # Computation of Gabor coefficient of dual frame window.
    dgrn1 = np.zeros((signal.shape[0], n_coeff), dtype=complex)
    k = np.arange(1, signal.shape[0] + 1)
    for n in range(n_coeff):
        index = modulo(k - n * m / q_oversample, signal.shape[0]).astype(int) - 1
        dgrn1[:, n] = np.fft.fft(signal * np.fft.fftshift(gam[index]), axis=0)
    dgr = dgrn1[np.arange(signal.shape[0], step=nb).astype(int), :]
    tfr = np.abs(dgr) ** 2
    return tfr, dgr, gam


# def original_main():
#     from tftb.generators import fmconst
#     import matplotlib.pyplot as plt
#     sig = np.r_[fmconst(128, 0.2)[0], fmconst(128, 0.4)[0]]
#     ts = np.linspace(0, 1, 256)
#     tfr = ShortTimeFourierTransform(sig, timestamps=ts)
#     tfr.run()
#     tfr.plot(show_tf=True, cmap=plt.cm.viridis)


def main():
    from tftb.generators import fmconst
    import matplotlib.pyplot as plt

    sig = np.r_[fmconst(128, 0.2)[0], fmconst(128, 0.4)[0]]
    ts = np.linspace(0, 1, 256)


    nr_samples = 256
    n_fbins = nr_samples
    nperseg =65
    noverlap = nperseg - 1
    window = hamming(nperseg)
    nfft = 256
    nyquist = nfft // 2
    print(f"nr_samples:  {nr_samples}")
    print(f"nperseg:  {nperseg}")
    print(f"noverlap:  {noverlap}")
    print(f"nfft:  {nfft}")
    stft = ShortTimeFourierTransform(sig, timestamps=None, n_fbins=n_fbins)
    tfr, ts, freqs = stft.run(
        nfft=nfft,
        nperseg=nperseg,
        noverlap=noverlap,
        return_onesided=False,
        window=window,
        scaling="psd")
    stft.plot(show_tf=True, cmap=plt.cm.viridis)


if __name__ == "__main__":
    main()
    # original_main()