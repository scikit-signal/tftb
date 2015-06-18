#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 2.6 of the tutorial.
Doppler Signal.
"""

from tftb.generators import doppler
from numpy import real
import matplotlib.pyplot as plt

fm, am, _ = doppler(256.0, 200.0, 4000.0 / 60.0, 10.0, 50.0)
signal = am * fm

plt.plot(real(signal))
plt.xlabel('Time')
plt.ylabel('Real part')
plt.title('Doppler')
plt.xlim(0, 256)
plt.grid()
plt.show()
