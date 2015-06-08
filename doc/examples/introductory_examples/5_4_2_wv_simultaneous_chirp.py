#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Wigner Ville distribution of two simultaneous chirps.
"""

from tftb.generators.api import fmlin, sigmerge
from tftb.processing.cohen import wigner_ville
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt

N = 64
sig = sigmerge(fmlin(N, 0, 0.4)[0], fmlin(N, 0.3, 0.5)[0], 1)
tfr = wigner_ville(sig)
tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0

t = np.arange(N)
f = np.linspace(0, 0.5, N)

fig, axImage = plt.subplots()
axImage.contour(t, f, tfr, 5)
axImage.grid(True)
axImage.set_title("Wigner Ville Distribution")
axImage.set_ylabel('Frequency')
axImage.set_xlabel('Time')

divider = make_axes_locatable(axImage)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(sig))
axTime.set_xlim(0, 64)
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axTime.grid(True)
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1][:32], np.arange(32))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
