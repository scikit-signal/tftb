#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example showing a Unterberger distribution of a hyperbolic group delay signal.
"""

from tftb.processing import unterberger
from tftb.generators import gdpower
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

sig = gdpower(128, -1)[0]
tfr, t, f = unterberger(sig, fmin=0.01, fmax=0.22, n_voices=172)
tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0
t, f = np.meshgrid(t, f)

fig, axContour = plt.subplots()
axContour.contour(t, f, tfr)
axContour.grid(True)
axContour.set_title("Unterberger distribution of hyperbolic GD signal.")
axContour.set_ylabel('Frequency')
axContour.set_xlabel('Time')

divider = make_axes_locatable(axContour)
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
plt.show()
