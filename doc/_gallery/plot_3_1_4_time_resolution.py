#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
=====================
Ideal time resolution
=====================

This example demonstrates that only the shortest possible window can provide
ideal resolution in time.

Figure 3.5 from the tutorial.
"""

import numpy as np
from tftb.processing.linear import ShortTimeFourierTransform
from tftb.generators import fmlin, amgauss
from matplotlib.pyplot import cm

x = np.real(amgauss(128) * fmlin(128)[0])
window = np.array([1])

nr_samples = 128
n_fbins = nr_samples
nperseg = 1
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
stft.plot(show_tf=True, cmap=cm.gray)
