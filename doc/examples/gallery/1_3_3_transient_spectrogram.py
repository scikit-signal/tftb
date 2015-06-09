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
from scipy.signal import hamming
from tftb.generators.api import amexpos, fmconst, sigmerge, noisecg
from tftb.processing.cohen import spectrogram

# Generate a noisy transient signal.
transsig = amexpos(64, kind='unilateral') * fmconst(64)[0]
signal = np.hstack((np.zeros((100,)), transsig, np.zeros((92,))))
signal = sigmerge(signal, noisecg(256), -5)

fwindow = hamming(65)

tfr, _, _ = spectrogram(signal, n_fbins=128, window=fwindow)
threshold = np.amax(tfr) * 0.1
tfr[tfr <= threshold] = 0.0
t = np.arange(256)
f = np.linspace(0, 0.5, 128)
t, f = np.meshgrid(t, f)

plt.contour(t, f, tfr)
plt.show()
