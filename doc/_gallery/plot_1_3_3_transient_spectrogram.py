#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
=======================================
Spectrogram of a Noisy Transient Signal
=======================================

This example demonstrates the simple use of a Spectrogram to localize a signal
in time and frequency. The transient signal appears at the normalized frequency
0.25 and between time points 125 and 160.

Figure 1.11 from the tutorial.
"""


import numpy as np
from scipy.signal.windows import hamming
from tftb.generators import amexpos, fmconst, sigmerge, noisecg
from tftb.processing.cohen import Spectrogram

# Generate a noisy transient signal.
transsig = amexpos(64, kind='unilateral') * fmconst(64)[0]
signal = np.hstack((np.zeros((100,)), transsig, np.zeros((92,))))
signal = sigmerge(signal, noisecg(256), -5)

fs = 1.0
nfft = 128
nperseg = 65
window = hamming(nperseg)
noverlap = nperseg - 1
# detrend = 'constant'
detrend = False
return_onesided = False
scaling = 'density'
mode = 'psd'

spec = Spectrogram(signal, n_fbins=nfft, fwindow=window)
spec.run(fs=1.0,
         window=window,
         nperseg=nperseg,
         noverlap=noverlap,
         nfft=nfft,
         # detrend=detrend,
         return_onesided=return_onesided,
         scaling=scaling,
         mode=mode
         )

spec.plot(kind="contour", threshold=0.1, show_tf=False)
