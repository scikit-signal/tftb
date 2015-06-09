#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example of the reassigned spectrogram.
"""

from tftb.generators.api import fmsin, fmhyp
from tftb.processing.api import ideal_tfr, reassigned_spectrogram, spectrogram
import numpy as np
import matplotlib.pyplot as plt

n_points = 128
sig1, if1 = fmsin(n_points, 0.15, 0.45, 100, 1, 0.4, -1)
sig2, if2 = fmhyp(n_points, [1, .5], [32, 0.05])
sig = sig1 + sig2
ideal, t, f = ideal_tfr(np.vstack((if1, if2)))
_, re_spec, _ = reassigned_spectrogram(sig)
spec, t3, f3 = spectrogram(sig)

# Ideal tfr
plt.subplot(221)
plt.contour(t, f, ideal, 1)
plt.grid(True)
plt.title("Ideal time-frequency distro")
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')

# Spectrogram
plt.subplot(222)
plt.contour(t3, f3[:64], spec[:64, :])
plt.grid(True)
plt.title("Spectrogram")
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')

# Reassigned Spectrogram
plt.subplot(212)
f = np.linspace(0, 0.5, 64)
plt.contour(np.arange(128), f, re_spec[:64, :])
plt.grid(True)
plt.title("Reassigned Spectrogram")
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')

plt.show()
