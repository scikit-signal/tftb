#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
==================================================
Linear Frequency and Gaussian Amplitude Modulation
==================================================

Generate a mono-component nonstationary signal with linear frequency
modulation and Gaussian amplitude modulation.

"""

from tftb.generators import fmlin, amgauss
from numpy import real
import matplotlib.pyplot as plt


fm, _ = fmlin(256)
am = amgauss(256)
signal = fm * am
plt.plot(real(signal))
plt.xlabel('Time')
plt.ylabel('Real part')
plt.title('Linear Frequency, Gaussian Amplitude')
plt.xlim(0, 256)
plt.grid()
plt.show()
