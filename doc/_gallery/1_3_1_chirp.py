#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 1.3.1 of the tutorial.
"""

from tftb.generators.api import fmlin, sigmerge, noisecg
from tftb.processing.cohen import wigner_ville
import matplotlib.pyplot as plt
import numpy as np

# Generate a chirp signal

n_points = 128
fmin, fmax = 0.0, 0.5

signal, _ = fmlin(n_points, fmin, fmax)
plt.plot(np.real(signal))
plt.title('Linear Frequency Modulation')
plt.ylabel('Real Part')
plt.xlabel('Time')
plt.grid()
plt.show()

# Plot the energy spectrum of the chirp

dsp1 = np.fft.fftshift(np.abs(np.fft.fft(signal)) ** 2)
plt.plot(np.arange(-64, 64, dtype=float) / 128.0, dsp1)
plt.title('Spectrum')
plt.ylabel('Squared modulus')
plt.xlabel('Normalized Frequency')
plt.grid()
plt.show()

# Wigner-Ville distribution of the chirp.

tfr = wigner_ville(signal)
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0
plt.contour(tfr, extent=[0, n_points, fmin, fmax])
plt.grid()
plt.title('Wigner-Ville distribution (threshold = 5%)')
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')
plt.show()

# Noisy chirp

noisy_signal = sigmerge(signal, noisecg(128), 0)
plt.plot(np.real(noisy_signal))
plt.title('Noisy chirp')
plt.ylabel('Real Part')
plt.xlabel('Time')
plt.grid()
plt.show()

# Enery spectrum of the noisy chirp.

dsp1 = np.fft.fftshift(np.abs(np.fft.fft(noisy_signal)) ** 2)
plt.plot(np.arange(-64, 64, dtype=float) / 128.0, dsp1)
plt.title('Spectrum of Noisy Chirp')
plt.ylabel('Squared modulus')
plt.xlabel('Normalized Frequency')
plt.grid()
plt.show()

# Wigner-Ville spectrum of noisy chirp.

tfr = wigner_ville(noisy_signal)
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0
plt.contour(tfr, extent=[0, n_points, fmin, fmax])
plt.grid()
plt.title('Wigner-Ville distribution of Noisy Signal (threshold = 5%)')
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')
plt.show()
