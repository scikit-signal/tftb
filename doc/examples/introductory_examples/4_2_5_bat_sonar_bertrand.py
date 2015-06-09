#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Bertrand distribution of Bat sonar
"""

from scipy.io import loadmat
from scipy.signal import hilbert
from tftb.processing.affine import bertrand
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

bat = loadmat("data/bat.mat")['bat']
N = 128
sig = hilbert(bat[np.arange(801, 800 + N * 7 + 1, step=7)])
tfr, t, f = bertrand(sig.ravel(), fmin=0.1, fmax=0.4, n_voices=N)

tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.01
tfr[tfr <= threshold] = 0.0

t, f = np.meshgrid(t, f)

fig, axImage = plt.subplots()
axImage.contour(t, f, tfr)
axImage.grid(True)
axImage.set_title("Bertrand distribution of bat sonar")
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
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2),
            np.arange(sig.shape[0]))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
