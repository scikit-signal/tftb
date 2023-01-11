#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
==================================================================
Spectrogram of Parallel Chirps with Short Gaussian Analysis Window
==================================================================

This example visualizes the spectrogram of two "parallel" chirps, using a
Gaussian window function that has a short length, relative to the length of a
signal. The two chirps can be made out, but interference can also be seen along
the frequency axis, since frequency resolution is compromised.

Figure 3.15 from the tutorial.
"""

from tftb.generators import fmlin
from tftb.processing.cohen import Spectrogram
import numpy as np
import matplotlib.pyplot as plt

sig = fmlin(128, 0, 0.4)[0] + fmlin(128, 0.1, 0.5)[0]
window = np.exp(np.log(0.005) * np.linspace(-1, 1, 23) ** 2)
spec = Spectrogram(sig, fwindow=window, n_fbins=128)
spec.run()
spec.plot(show_tf=True, cmap=plt.cm.gray)
