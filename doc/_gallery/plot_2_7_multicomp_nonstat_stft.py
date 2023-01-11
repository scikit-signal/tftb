#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
======================================================================
Short time Fourier transform of a multi-component nonstationary signal
======================================================================

Compute and visualize the `STFT <https://en.wikipedia.org/wiki/Short-time_Fourier_transform>`_ of a multi component nonstationary signal.

Figure 2.11 from the tutorial.
"""

from tftb.generators import fmlin
# from tftb.processing.linear import ShortTimeFourierTransform
from tftb.processing.linear import ShortTimeFourierTransform as ShortTimeFourierTransform
import matplotlib.pyplot as plt
from scipy.signal.windows import hamming
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


def generate_signal(nr_samples):
    x1, _ = fmlin(nr_samples, 0, 0.2)
    x2, _ = fmlin(nr_samples, 0.3, 0.5)
    return x1 + x2


def plot_stft_contours(axis, ts, freqs, tfr):
    print("In plot_stft_contours")
    print(f"freqs.shape:  {freqs.shape}")
    print(freqs)
    t_mesh, f_mesh = np.meshgrid(ts, freqs)
    # axis.contour(T[:32, :], F[:32], data[:32, :], 5)
    nyquist = freqs.shape[0] // 2
    print(f"nyquist: {nyquist}")
    axis.contour(t_mesh[:nyquist, :], f_mesh[:nyquist], tfr[:nyquist, :], 5)
    axis.grid(True)
    axis.set_title('Squared modulus of STFT')
    axis.set_ylabel('Frequency')
    axis.yaxis.set_label_position("right")
    axis.set_xlabel('Time')
    return axis


def plot_signal(axis, x, nr_samples):
    axis.plot(np.real(x))
    axis.set_xticklabels([])
    axis.set_xlim(0, nr_samples)
    axis.set_ylabel('Real part')
    axis.set_title('Signal in time')
    axis.grid(True)


def plot_spectrum(axis, squared_modulus_spectrum, freqs, nyquist):
    print("In plot_spectrum")
    print(f"freqs.shape:  {freqs.shape}")
    print(freqs)
    print(f"squared_modulus_spectrum.shape:  {squared_modulus_spectrum.shape}")
    print(f"max: {np.max(freqs)}")
    print(f"min: {np.min(freqs)}")
    len_freqs = squared_modulus_spectrum.shape[0] // 2
    # len_freqs = freqs.shape[0] // 2
    print(f"len_freqs: {len_freqs}")
    print(f"nyquist: {nyquist}")
    # axis.plot(squared_modulus_spectrum[::-1][:nyquist], freqs[:nyquist])
    f = np.linspace(0.0, 0.5, len_freqs, endpoint=True)
    axis.plot(squared_modulus_spectrum[len_freqs:], f)
    axis.set_yticklabels([])
    axis.set_xticklabels([])
    axis.invert_xaxis()
    axis.set_ylabel('Spectrum')
    axis.grid(True)


def plot(nr_samples, x, tfr, ts, freqs, squared_modulus_spectrum, nyquist):
    fig, axScatter = plt.subplots(figsize=(10, 8))

    plot_stft_contours(axScatter, ts, freqs, tfr)

    divider = make_axes_locatable(axScatter)
    axTime = divider.append_axes("top", 1.2, pad=0.5)
    axFreq = divider.append_axes("left", 1.2, pad=0.5)

    plot_signal(axTime, x, nr_samples)
    plot_spectrum(axFreq, squared_modulus_spectrum, freqs, nyquist)
    plt.show()


def main():
    nr_samples = 128
    x = generate_signal(nr_samples)
    n_fbins = nr_samples
    nperseg = 33
    noverlap = nperseg - 1
    window = hamming(nperseg)
    nfft = 128
    nyquist = nfft // 2
    print(f"nr_samples:  {nr_samples}")
    print(f"nperseg:  {nperseg}")
    print(f"noverlap:  {noverlap}")
    print(f"nfft:  {nfft}")
    stft = ShortTimeFourierTransform(x, timestamps=None, n_fbins=n_fbins)
    tfr, ts, freqs = stft.run(
        nfft=nfft,
        nperseg=nperseg,
        noverlap=noverlap,
        return_onesided=False,
        window=window,
        scaling="psd")

    print(f"tfr shape:  {tfr.shape}")
    # data, _, _ = ShortTimeFourierTransform(x, timestamps=None, n_fbins=n_fbins,
    #                                       fwindow=window).run()
    # data = data[:64, :]
    print(f"tfr shape modified:  {tfr.shape}")
    threshold = np.amax(np.abs(tfr)) * 0.05
    print(f"threshold:  {threshold}")

    tfr[np.abs(tfr) <= threshold] = 0.0 + 1j * 0.0
    tfr = np.abs(tfr) ** 2
    print(f"len(x): {len(x)}")
    squared_modulus_spectrum = abs(np.fft.fftshift(np.fft.fft(x))) ** 2
    print(f"In main. squared_modulus_spectrum:  {squared_modulus_spectrum.shape}")
    plot(nr_samples, x, tfr, ts, freqs, squared_modulus_spectrum, nyquist)


if __name__ == "__main__":
    main()
