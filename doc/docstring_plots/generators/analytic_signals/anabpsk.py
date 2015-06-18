#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import anabpsk
import matplotlib.pyplot as plt
import numpy as np

x, am = anabpsk(300, 30, 0.1)
plt.subplot(211), plt.plot(np.real(x))
plt.grid()
plt.title('Analytic BPSK signal')
plt.subplot(212), plt.plot(am)
plt.grid()
plt.title('Amplitude Modulation')
plt.show()
