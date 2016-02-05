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
stft = ShortTimeFourierTransform(x, n_fbins=128, fwindow=window)
stft.run()
stft.plot(show_tf=True, cmap=plt.cm.gray)
