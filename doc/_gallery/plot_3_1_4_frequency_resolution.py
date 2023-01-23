#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
================================================
Time-frequency Resolution: Long Analysis Window
================================================

This example shows the effect of an analysis window which is long in time on
the time-frequency resolution. Specifically, longer windows have good frequency
resolutions but poor time resolutions.

Figure 3.6 from the tutorial.
"""

import numpy as np
from tftb.processing.linear import ShortTimeFourierTransform
from tftb.generators import fmlin, amgauss
import matplotlib.pyplot as plt

x = np.real(amgauss(128) * fmlin(128)[0])
window = np.ones((128,))
nr_samples = 128
n_fbins = nr_samples
nperseg = 128
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
