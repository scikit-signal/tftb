#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 1.3.3 of the tutorial.
"""


# dsp=fftshift(abs(fft(sign)).^2);
# plot((-128:127)/256,dsp);

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators.api import amexpos, fmconst, sigmerge, noisecg
from tftb.processing.reassigned import spectrogram
from scipy.signal import hamming

# Generate a noisy transient signal.
transsig = amexpos(64, kind='unilateral') * fmconst(64)[0]
signal = np.hstack((np.zeros((100,)), transsig, np.zeros((92,))))
signal = sigmerge(signal, noisecg(256), -5)
plt.plot(np.real(signal))
plt.grid()
plt.title('Noisy Transient Signal')
plt.xlabel('Time')
plt.xlim((0, 256))
plt.ylim((np.real(signal).max(), np.real(signal.min())))
plt.show()

# Energy spectrum of the signal
dsp = np.fft.fftshift(np.abs(np.fft.fft(signal)) ** 2)
plt.plot(np.arange(-128, 128, dtype=float) / 256, dsp)
plt.title('Energy spectrum of noisy transient signal')
plt.xlabel('Normalized frequency')
plt.grid()
plt.xlim(-0.5, 0.5)
plt.show()

# Get reassinged spectrogram of the signal.
window = hamming(65)
tfr, _, _ = spectrogram(signal, n_fbins=128, window=window)
threshold = np.amax(tfr) * 0.1
tfr[tfr <= threshold] = 0.0
f = np.linspace(0, 0.5, 128)
t = np.arange(256)
plt.contour(t, f, tfr, 5)
plt.grid()
plt.show()
