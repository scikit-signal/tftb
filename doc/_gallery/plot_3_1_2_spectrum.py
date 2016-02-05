#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
==================================
Energy Spectrum of an Audio Signal
==================================

Figure 3.3 from the tutorial.
"""

from os.path import dirname, abspath, join
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = join(abspath(dirname("__file__")), "data", "gabor.mat")
signal = loadmat(DATA_PATH)['gabor'].ravel()
time = np.arange(338)
dsp = np.fft.fftshift(np.abs(np.fft.fft(signal)) ** 2)
freq = np.arange(-169, 169, dtype=float) / 338 * 1000

plt.subplot(211)
plt.plot(time, signal)
plt.grid(True)
plt.xlim(0, time.max())
plt.xlabel('Time (ms)')

plt.subplot(212)
plt.plot(dsp)
plt.grid(True)
plt.title('Spectrum')
plt.xlabel('Frequency (Hz)')

plt.subplots_adjust(hspace=0.5)

plt.show()
