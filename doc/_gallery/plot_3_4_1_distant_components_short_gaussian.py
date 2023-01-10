#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
====================================================
Distant Chirps with a Short Gaussian Analysis Window
====================================================

Figure 3.17 from the tutorial.
"""

from tftb.generators import fmlin
from tftb.processing.cohen import Spectrogram
import numpy as np
import matplotlib.pyplot as plt

sig = fmlin(128, 0, 0.3)[0] + fmlin(128, 0.2, 0.5)[0]
window = np.exp(np.log(0.005) * np.linspace(-1, 1, 23) ** 2)

nperseg = 23
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
