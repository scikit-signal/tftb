#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 2.7 of the tutorial.
Short time Fourier transform of a multi-component nonstationary signal.
"""

from tftb.generators.api import fmlin
from tftb.processing.stft import stft
import matplotlib.pyplot as plt
from scipy.signal import hamming
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

N = 128
x1, _ = fmlin(N, 0, 0.2)
x2, _ = fmlin(N, 0.3, 0.5)
x = x1 + x2

n_fbins = 128
window = hamming(33)
tfr, _, _ = stft(x, None, n_fbins, window)
threshold = np.amax(np.abs(tfr)) * 0.05
tfr[np.abs(tfr) <= threshold] = 0.0 + 1j * 0.0
tfr = np.abs(tfr) ** 2
t = np.arange(tfr.shape[1])
f = np.linspace(0, 0.5, tfr.shape[0])

fig, axScatter = plt.subplots()
axScatter.contour(tfr, 5)
axScatter.grid(True)
axScatter.set_title('Squared modulus of STFT')
axScatter.set_ylabel('Frequency')
axScatter.set_xlabel('Time')
divider = make_axes_locatable(axScatter)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(x))
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axTime.grid(True)
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(x))) ** 2)[::-1], np.arange(x.shape[0]))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
