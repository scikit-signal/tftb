#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example showing a Morlet scalogram of two atoms.
"""

from tftb.processing.api import scalogram
from tftb.generators import atoms
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

# FIXME: squared modulus of scalograms doesn't seem to be necessary.

sig = atoms(128, np.array([[38, 0.1, 32, 1], [96, 0.35, 32, 1]]))
tfr, t, f, _ = scalogram(sig)
t, f = np.meshgrid(t, f)

fig, axContour = plt.subplots()
axContour.contour(t, f, tfr)
axContour.grid(True)
axContour.set_title("Morlet scalogram of complex sinusoids")
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
