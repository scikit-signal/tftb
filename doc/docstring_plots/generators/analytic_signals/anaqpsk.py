#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import anaqpsk
import matplotlib.pyplot as plt
import numpy as np

x, am = anaqpsk(512, 64.0, 0.05)
plt.subplot(211), plt.plot(np.real(x))
plt.xlim(0, 512)
plt.grid()
plt.title('Analytic QPSK signal')
plt.subplot(212), plt.plot(am)
plt.xlim(0, 512)
plt.grid()
plt.title('Phase')
plt.show()
