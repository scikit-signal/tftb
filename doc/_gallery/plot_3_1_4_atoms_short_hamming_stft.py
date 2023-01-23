#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
================================================
Time-frequency Resolution: Short Analysis Window
================================================

This example shows the effect of an analysis window which is short in time on
the time-frequency resolution. Specifically, smaller windows have good time
resolutions but poor frequency resolutions.

Figure 3.8 from the tutorial.
"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators import atoms
from scipy.signal import hamming
from tftb.processing.linear import ShortTimeFourierTransform


coords = np.array([[45, .25, 32, 1], [85, .25, 32, 1]])
sig = atoms(128, coords)
x = np.real(sig)

nperseg = 17
window = hamming(nperseg)

n_fbins = 128
noverlap = nperseg - 1
nfft = 128
stft = ShortTimeFourierTransform(x, timestamps=None, n_fbins=n_fbins)
tfr, ts, freqs = stft.run(
    nfft=nfft,
    nperseg=nperseg,
    noverlap=noverlap,
    return_onesided=False,
    window=window,
    scaling="psd")
stft.plot(show_tf=True, cmap=plt.cm.gray)
