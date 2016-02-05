#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
==============================================================
Narrow Band Ambiguity Function of Chirps with Different Slopes
==============================================================

This example demonstrates the narrow band ambiguity function (AF) of a signal
composed of two chirps with Gaussian amplitude modulation but havind linear
frequency modulations with different slopes. Note that the AF interference
terms are located away from the origin.

Figure 4.13 from the tutorial.
"""

from tftb.generators import fmlin, amgauss
from tftb.processing.ambiguity import narrow_band
import numpy as np
import matplotlib.pyplot as plt

n_points = 64
sig1 = fmlin(n_points, 0.2, 0.5)[0] * amgauss(n_points)
sig2 = fmlin(n_points, 0.3, 0)[0] * amgauss(n_points)
sig = np.hstack((sig1, sig2))

tfr, x, y = narrow_band(sig)
plt.contour(2 * x, y, np.abs(tfr) ** 2, 16)
plt.title('Narrow Band ambiguity function')
plt.xlabel('Delay')
plt.ylabel('Doppler')
plt.grid(True)
plt.show()
