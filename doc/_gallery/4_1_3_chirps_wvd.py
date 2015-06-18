#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Examples from section 4.1.3 of the tutorial.
"""

from tftb.generators import fmlin, amgauss
from tftb.processing.api import wigner_ville
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

n_points = 64
sig1 = fmlin(n_points, 0.2, 0.5)[0] * amgauss(n_points)
sig2 = fmlin(n_points, 0.3, 0)[0] * amgauss(n_points)
sig = np.hstack((sig1, sig2))

tfr = wigner_ville(sig)
x = np.arange(tfr.shape[0])
y = np.linspace(0, 0.5, tfr.shape[0])

threshold = np.amax(np.abs(tfr)) * 0.05
tfr[np.abs(tfr) <= threshold] = 0.0

T, F = np.meshgrid(x, y)

fig, axContour = plt.subplots(figsize=(10, 8))
axContour.contour(T, F, tfr)
axContour.grid(True)
axContour.set_title("Wigner Ville distribution")
axContour.set_ylabel('Frequency')
axContour.yaxis.set_label_position('right')
axContour.set_xlabel('Time')

divider = make_axes_locatable(axContour)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(sig))
axTime.set_xticklabels([])
axTime.set_xlim(0, 128)
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axTime.grid(True)
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1][:64],
            np.arange(sig.shape[0] / 2))
axFreq.set_ylim(0, sig.shape[0] / 2 - 1)
axFreq.set_yticklabels([])
axFreq.set_xticklabels([])
axFreq.grid(True)
axFreq.set_ylabel('Spectrum')
axFreq.invert_xaxis()
axFreq.grid(True)
plt.show()
