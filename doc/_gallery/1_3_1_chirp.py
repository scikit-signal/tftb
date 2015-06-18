#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 1.3.1 of the tutorial.
"""

from tftb.generators import fmlin
import matplotlib.pyplot as plt
import numpy as np

# Generate a chirp signal

n_points = 128
fmin, fmax = 0.0, 0.5

signal, _ = fmlin(n_points, fmin, fmax)
plt.plot(np.real(signal))
plt.xlim(0, n_points)
plt.title('Linear Frequency Modulation')
plt.ylabel('Real Part')
plt.xlabel('Time')
plt.grid()
plt.show()
