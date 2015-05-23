#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 2.2.2 of the tutorial.
"""

from tftb.generators.api import amgauss
from tftb.processing.api import loctime, locfreq
import numpy as np
import matplotlib.pyplot as plt

# generate signal
signal = amgauss(256)
plt.subplot(211), plt.plot(np.real(signal))
plt.xlim(0, 256)
plt.xlabel('Time')
plt.ylabel('Real part')
plt.title('Signal')
plt.grid()
fsig = np.fft.fftshift(np.abs(np.fft.fft(signal)) ** 2)
plt.subplot(212), plt.plot(np.linspace(0, 0.5, 256), fsig)
plt.xlabel('Normalized frequency')
plt.ylabel('Squared modulus')
plt.title('Spectrum')
plt.grid()
plt.show()

tm, T = loctime(signal)
print "Time Center: {}".format(tm)
print "Time Duration: {}".format(T)
fm, B = locfreq(signal)
print "Frequency Center: {}".format(fm)
print "Frequency Spreading: {}".format(B)
print "Time-bandwidth product: {}".format(T * B)
