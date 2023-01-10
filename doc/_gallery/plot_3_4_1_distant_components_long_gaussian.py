#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
===================================================
Distant Chirps with a Long Gaussian Analysis Window
===================================================

This example visualizes a spectrogram of two chirp signals which are well
separated in frequency ranges. A longer Gaussian analysis window suffices here
to see the separation of frequencies, since the variation in frequencies is
relatively slow.

Figure 3.18 from the tutorial.
"""

from tftb.generators import fmlin
from tftb.processing.cohen import Spectrogram
import numpy as np
import matplotlib.pyplot as plt

sig = fmlin(128, 0, 0.3)[0] + fmlin(128, 0.2, 0.5)[0]
window = np.exp(np.log(0.005) * np.linspace(-1, 1, 63) ** 2)

nperseg = 63
fs = 1.0
nfft = 128

window = np.exp(np.log(0.005) * np.linspace(-1, 1, nperseg) ** 2)
noverlap = nperseg - 1
# detrend = 'constant'
detrend = "constant"
return_onesided = False
scaling = 'spectrum'
mode = 'psd'

spec = Spectrogram(sig, n_fbins=nfft, fwindow=window)
spec.run(fs=1.0,
         window=window,
         nperseg=nperseg,
         noverlap=noverlap,
         nfft=nfft,
         detrend=detrend,
         return_onesided=return_onesided,
         scaling=scaling,
         mode=mode
         )

spec.plot(show_tf=True, cmap=plt.cm.gray)
