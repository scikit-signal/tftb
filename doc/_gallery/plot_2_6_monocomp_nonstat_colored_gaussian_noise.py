#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
=========================
Noisy Monocomponent Chirp
=========================

This example demonstrates the construction of a monocomponent signal with
linear frequency modulation and colored Gaussian noise.

Figure 2.9 from the tutorial.
"""

from tftb.generators import fmlin, amgauss, noisecg, sigmerge
from numpy import real
import matplotlib.pyplot as plt

fm, _ = fmlin(256)
am = amgauss(256)
signal = fm * am

noise = noisecg(256, .8)
sign = sigmerge(signal, noise, -10)

plt.plot(real(sign))
plt.xlabel('Time')
plt.ylabel('Real part')
plt.title('Gaussian transient signal embedded in -10 dB colored Gaussian noise')
plt.xlim(0, 256)
plt.grid()
plt.show()
