#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example showing use of Hough transform on a Wigner-Ville distribution.
"""

import numpy as np
from tftb.generators import noisecg, sigmerge, fmlin
from tftb.processing.cohen import wigner_ville
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

N = 64
sig = sigmerge(fmlin(N, 0, 0.3)[0], noisecg(N), 1)
tfr = wigner_ville(sig)
tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0

time = np.arange(N)
f = np.linspace(0, 0.5, N)

fig, axContour = plt.subplots(figsize=(10, 8))
axContour.contour(time, f, tfr, 5)
axContour.grid(True)
axContour.set_title("Wigner Ville Distribution")
axContour.yaxis.set_label_position('right')
axContour.set_ylabel('Frequency')
axContour.set_xlabel('Time')

divider = make_axes_locatable(axContour)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(sig))
axTime.set_xticklabels([])
axTime.set_xlim(0, 64)
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axFreq.set_ylim(0, sig.shape[0] / 2 - 1)
axFreq.set_yticklabels([])
axFreq.set_xticklabels([])
axTime.grid(True)
axFreq.invert_xaxis()
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1][:32], np.arange(32))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
