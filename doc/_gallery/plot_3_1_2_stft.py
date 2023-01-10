#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
=======================
STFT of an Audio Signal
=======================

Figure 3.4 from the tutorial.
"""

from os.path import dirname, abspath, join
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = join(abspath(dirname("__file__")), "data", "gabor.mat")
signal = loadmat(DATA_PATH)['gabor'].ravel()
tfr = loadmat(DATA_PATH)['tfr']
time = np.arange(338)
freq = np.arange(128, dtype=float) / 256.0 * 1000

plt.contour(time, freq, tfr)
plt.grid(True)
plt.xlabel('Time  [ms]')
plt.ylabel('Frequency  [Hz]')
plt.title('Squared modulus of the STFT of the word GABOR')
plt.show()
