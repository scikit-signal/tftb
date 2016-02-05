#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
===========================
Linear Frequency Modulation
===========================

This example shows how PyTFTB is used to generate a signal with linear
frequency modulation. Such a signal is also called a `chirp
<https://en.wikipedia.org/wiki/Chirp>`_.

Figure 1.1 from the tutorial.
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
