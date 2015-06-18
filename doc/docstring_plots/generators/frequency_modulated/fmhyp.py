#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.
"""

"""

from tftb.generators import fmhyp
import numpy as np
import matplotlib.pyplot as plt

signal, iflaw = fmhyp(128, (1, 0.5), (32, 0.1))
plt.subplot(211), plt.plot(np.real(signal))
plt.xlim(0, 128)
plt.grid()
plt.title('Hyperbolic Frequency Modulation')
plt.subplot(212), plt.plot(iflaw)
plt.xlim(0, 128)
plt.grid()
plt.title('Instantaneous Frequency')
plt.show()
