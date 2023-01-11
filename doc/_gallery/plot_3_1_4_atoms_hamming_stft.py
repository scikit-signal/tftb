#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
============================================================
STFT of Gaussian Wave Packets with a Hamming Analysis Window
============================================================

This example demonstrates the construction of a signal containing two transient
components, having the same Gaussian amplitude modulation and the same
frequency, but different time centers. It also shows the effect of a Hamming
window function when used with th STFT.

Figure 3.7 from the tutorial.
"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators import atoms
from scipy.signal.windows import hamming
from tftb.processing.linear import ShortTimeFourierTransform as ShortTimeFourierTransform

coords = np.array([[45, .25, 32, 1], [85, .25, 32, 1]])
sig = atoms(128, coords)
x = np.real(sig)
window = hamming(65)
stft = ShortTimeFourierTransform(sig, n_fbins=128, fwindow=window)

nfft = 128
nperseg = 65
noverlap = nperseg - 1
window = hamming(nperseg)
n_fbins = nfft
tfr, ts, freqs = stft.run(
    nfft=nfft,
    nperseg=nperseg,
    noverlap=noverlap,
    return_onesided=False,
    window=window,
    scaling="psd")
stft.plot(show_tf=True, cmap=plt.cm.gray)
stft.plot(show_tf=True, cmap=plt.cm.gray)
