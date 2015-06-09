#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators.api import anaask
import matplotlib.pyplot as plt
import numpy as np

x, am = anaask(512, 64, 0.05)
plt.subplot(211), plt.plot(np.real(x))
plt.grid()
plt.title('Analytic ASK signal')
plt.subplot(212), plt.plot(am)
plt.grid()
plt.title('Amplitude Modulation')
plt.show()
