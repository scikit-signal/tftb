#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 4.1.1 of the tutorial.
"""

import numpy as np
from tftb.generators.api import doppler
from tftb.processing.api import wigner_ville
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fm, am, iflaw = doppler(256, 50.0, 13.0, 10.0, 200.0)
sig = am * fm
tfr = wigner_ville(sig)
threshold = (np.abs(tfr) ** 2) * 0.05
tfr[np.abs(tfr) ** 2 <= threshold] = 0.0


fig, axImage = plt.subplots()
axImage.contour(np.flipud(tfr[:128, :]), extent=[0, 128, 0, 0.5])
axImage.grid(True)
axImage.set_title("TFRWV")
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
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1][:128],
            np.arange(sig.shape[0] / 2))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
