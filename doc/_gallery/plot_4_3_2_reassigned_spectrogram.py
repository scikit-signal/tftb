#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
========================================================
Comparison of a Spectrogram and a Reassigned Spectrogram
========================================================

This example compares the spectrogram and the reassigned spectrogram of a
hybrid signal (containing sinusoidal, constant and linear frequency
modulations), against its ideal time-frequency characteristics.

Figure 4.34 from the tutorial.
"""
from scipy.signal.windows import hamming

from tftb.processing.cohen import Spectrogram

from tftb.generators import fmsin, fmhyp
from tftb.processing import ideal_tfr, reassigned_spectrogram
import numpy as np
import matplotlib.pyplot as plt

n_points = 128
sig1, if1 = fmsin(n_points, 0.15, 0.45, 100, 1, 0.4, -1)
sig2, if2 = fmhyp(n_points, [1, .5], [32, 0.05])
sig = sig1 + sig2
ideal, t, f = ideal_tfr(np.vstack((if1, if2)))
_, re_spec, _ = reassigned_spectrogram(sig)

fs = 1.0
nfft = 128
nperseg = 33
window = hamming(nperseg)
noverlap = nperseg - 1
# detrend = 'constant'
detrend = "constant"
return_onesided = False
scaling = 'spectrum'
mode = 'psd'

spec = Spectrogram(sig, n_fbins=nfft, fwindow=window)
spec, t3, f3 = spec.run(
    fs=1.0,
    window=window,
    nperseg=nperseg,
    noverlap=noverlap,
    nfft=nfft,
    detrend=detrend,
    return_onesided=return_onesided,
    scaling=scaling,
    mode=mode
)

# Ideal data
plt.subplot(221)
plt.contour(t, f, ideal, 1)
plt.grid(True)
plt.gca().set_xticklabels([])
plt.title("Ideal time-frequency distro")
plt.ylabel('Normalized Frequency')

# Spectrogram
plt.subplot(222)
plt.contour(t3, f3[:64], spec[:64, :])
plt.grid(True)
plt.gca().set_xticklabels([])
plt.title("Spectrogram")

# Reassigned Spectrogram
plt.subplot(212)
f = np.linspace(0, 0.5, 64)
plt.contour(np.arange(128), f, re_spec[:64, :])
plt.grid(True)
plt.title("Reassigned Spectrogram")
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')

plt.show()
