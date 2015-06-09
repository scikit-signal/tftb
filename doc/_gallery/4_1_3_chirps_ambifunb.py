#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Examples from section 4.1.3 of the tutorial.
"""

from tftb.generators.api import fmlin, amgauss
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
