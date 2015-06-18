#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 4.1.2 of the tutorials.
"""

from tftb.generators import fmconst, amgauss
from tftb.processing import pseudo_wigner_ville
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

sig = fmconst(128, 0.15)[0] + amgauss(128) * fmconst(128, 0.4)[0]
tfr = pseudo_wigner_ville(sig)
threshold = np.abs(tfr) * 0.05
tfr[np.abs(tfr) <= threshold] = 0.0

t = np.arange(128)
f = np.linspace(0, 0.5, 128)
t, f = np.meshgrid(t, f)
fig, axImage = plt.subplots(figsize=(10, 8))
axImage.contour(t, f, np.abs(tfr), extent=[0, 128, 0, 0.5])
axImage.grid(True)
axImage.set_title("Pseudo WV distribution.")
axImage.set_ylabel('Frequency')
axImage.yaxis.set_label_position('right')
axImage.set_xlabel('Time')

divider = make_axes_locatable(axImage)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(sig))
axTime.set_xticklabels([])
axTime.set_xlim(0, 128)
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axTime.grid(True)
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1][:64],
        np.arange(sig.shape[0] / 2)[::-1])
axFreq.set_ylim(0, sig.shape[0] / 2 - 1)
axFreq.set_yticklabels([])
axFreq.set_xticklabels([])
axFreq.grid(True)
axFreq.set_ylabel('Spectrum')
axFreq.invert_xaxis()
axFreq.grid(True)
plt.show()
