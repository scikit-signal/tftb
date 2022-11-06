#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
===============================================
Time and Frequency Localization Characteristics
===============================================

Generate a signal that has localized characteristics in both time and frequency
and compute the following estimates:

    * time center
    * time duration
    * frequency center
    * frequency spreading

Example 2.1 from the tutorial.
"""

from tftb.generators import fmlin, amgauss
from tftb.processing import loctime, locfreq
import numpy as np
import matplotlib.pyplot as plt

# generate signal
signal = fmlin(256)[0] * amgauss(256)
plt.subplot(211), plt.plot(np.real(signal))
plt.xlim(0, 256)
plt.xlabel('Time')
plt.ylabel('Real part')
plt.title('Signal')
plt.grid()
fsig = np.fft.fftshift(np.abs(np.fft.fft(signal)) ** 2)
plt.subplot(212), plt.plot(np.linspace(-0.5, 0.5, 256), fsig)
plt.xlabel('Normalized frequency')
plt.ylabel('Squared modulus')
plt.title('Spectrum')
plt.grid()
plt.subplots_adjust(hspace=0.5)
plt.show()


tm, T = loctime(signal)
print("Time Center: {}".format(tm))
print("Time Duration: {}".format(T))
num, B = locfreq(signal)
print("Frequency Center: {}".format(num))
print("Frequency Spreading: {}".format(B))
