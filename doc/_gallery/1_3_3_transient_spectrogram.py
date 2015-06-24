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
from scipy.signal import hamming
from tftb.generators import amexpos, fmconst, sigmerge, noisecg
from tftb.processing.cohen import Spectrogram

# Generate a noisy transient signal.
transsig = amexpos(64, kind='unilateral') * fmconst(64)[0]
signal = np.hstack((np.zeros((100,)), transsig, np.zeros((92,))))
signal = sigmerge(signal, noisecg(256), -5)

fwindow = hamming(65)

spec = Spectrogram(signal, n_fbins=128, fwindow=fwindow)
spec.run()
spec.plot(kind='contour')
