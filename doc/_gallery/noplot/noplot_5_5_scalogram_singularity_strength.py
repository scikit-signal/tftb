#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
===========================================
Morlet Scalogram of a Lipschitz Singularity
===========================================

The time localization of the Lipschitz function can be seen at smaller scales.

Figure 5.8 from the tutorial.
"""

from tftb.processing import scalogram
from tftb.generators import anasing
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt


sig = anasing(64, 32, -0.5)

tfr, t, f, _ = scalogram(sig, waveparams=4, fmin=0.01, fmax=0.5, n_voices=256)

t, f = np.meshgrid(t, f)

fig, axContour = plt.subplots()
axContour.contour(t, f, tfr, 10)
axContour.grid(True)
axContour.set_title("Morlet Scalogram of Lipschitz singularity")
axContour.set_ylabel('Frequency')
axContour.set_xlabel('Time')

divider = make_axes_locatable(axContour)
axTime = divider.append_axes("top", 1.2, pad=0.5)
axFreq = divider.append_axes("left", 1.2, pad=0.5)
axTime.plot(np.real(sig))
axTime.set_xlim(0, 64)
axTime.set_ylabel('Real part')
axTime.set_title('Signal in time')
axTime.grid(True)
axFreq.plot((abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1],
            np.arange(sig.shape[0]))
axFreq.set_ylabel('Spectrum')
axFreq.grid(True)
plt.show()
plt.show()
