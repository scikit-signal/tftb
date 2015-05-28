#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Examples from section 3.4.1 of the tutorial.
"""

from tftb.generators.api import fmlin
from tftb.processing.cohen import spectrogram
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

sig = fmlin(128, 0, 0.3)[0] + fmlin(128, 0.2, 0.5)[0]
window = np.exp(np.log(0.005) * np.linspace(-1, 1, 23) ** 2)
tfr, t, f = spectrogram(sig, window=window, n_fbins=128)
threshold = np.amax(np.abs(tfr)) * 0.05
tfr[np.abs(tfr) <= threshold] = 0.0

fig, axImage = plt.subplots()
axImage.imshow(np.flipud(tfr[:64, :]), extent=[0, 128, 0, 0.5], aspect='auto',
               cmap=plt.cm.gray)
axImage.grid(True)
axImage.set_title("Spectrogram of distant components - short Gaussian window")
axImage.set_ylabel('Frequency')
axImage.set_xlabel('Time')

divider = make_axes_locatable(axImage)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(sig))
axTime.set_xlim(0, 128)
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axTime.grid(True)
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1][:64],
            np.arange(sig.shape[0] / 2))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
