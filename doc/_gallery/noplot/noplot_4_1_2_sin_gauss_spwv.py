#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 4.1.2 of the tutorials.

Figure 4.10 from the tutorial.
"""

from tftb.generators import fmconst, amgauss
from tftb.processing import smoothed_pseudo_wigner_ville
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import kaiser, hamming

twindow = kaiser(13, 3 * np.pi)
fwindow = kaiser(33, 3 * np.pi)
twindow = hamming(13)
fwindow = hamming(33)

sig = fmconst(128, 0.15)[0] + amgauss(128) * fmconst(128, 0.4)[0]
tfr = smoothed_pseudo_wigner_ville(sig, twindow=twindow, fwindow=fwindow,
                                   freq_bins=128)
threshold = np.abs(tfr) * 0.05
tfr[np.abs(tfr) <= threshold] = 0.0

fig, axImage = plt.subplots()
axImage.contour(np.abs(tfr), extent=[0, 128, 0, 0.5], levels=list(range(5)))
axImage.grid(True)
axImage.set_title("Wigner Ville distribution")
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
spectrum = abs(np.fft.fftshift(np.fft.fft(sig)) ** 2)[64:]
axFreq.plot(spectrum, np.arange(spectrum.shape[0]))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
