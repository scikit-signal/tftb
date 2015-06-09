#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Examples from section 3.1.4 of the tutorial.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tftb.generators.api import atoms
from scipy.signal import hamming
from tftb.processing.linear import stft

coords = np.array([[45, .25, 32, 1], [85, .25, 32, 1]])
sig = atoms(128, coords)
x = np.real(sig)
window = hamming(65)
tfr, _, _ = stft(sig, n_fbins=128, window=window)
threshold = np.abs(tfr) * 0.05
tfr[np.abs(tfr) <= threshold] = 0.0 * 1j * 0.0

fig, axImage = plt.subplots()
mappable = axImage.imshow(np.abs(tfr[:64, :]) ** 2, cmap=plt.cm.gray)
fig.colorbar(mappable, ax=axImage, shrink=0.75)
axImage.grid(True)
axImage.set_title("STFT squared modulus")
axImage.set_ylabel('Frequency')
axImage.set_xlabel('Time')

divider = make_axes_locatable(axImage)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(x))
axTime.set_xlim(0, 128)
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axTime.grid(True)
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(x))) ** 2)[::-1][:64],
            np.arange(x.shape[0] / 2))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
