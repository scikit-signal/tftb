#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 2.6 of the tutorial.
Embedding a mono-component nonstationary signal with linear frequency
modulation and Gaussian amplitude modulation into Gaussian colored noise.
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
