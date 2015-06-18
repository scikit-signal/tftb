#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 2.6 of the tutorial.
Generating a mono-component nonstationary signal with constant frequency
modulation and one-sided exponential amplitude modulation.
"""

from tftb.generators import fmconst, amexpos
import matplotlib.pyplot as plt
from numpy import real

fm, _ = fmconst(256, 0.2)
am = amexpos(256, 100, kind='unilateral')
signal = am * fm

plt.plot(real(signal))
plt.xlabel('Time')
plt.ylabel('Real part')
plt.title('Constant Frequency, One-sided Exponential Amplitude')
plt.xlim(0, 256)
plt.grid()
plt.show()
