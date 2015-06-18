#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import anafsk
import matplotlib.pyplot as plt
import numpy as np

x, am = anafsk(512, 54.0, 5.0)
plt.subplot(211), plt.plot(np.real(x))
plt.xlim(0, 512)
plt.grid()
plt.title('Analytic FSK signal')
plt.subplot(212), plt.plot(am)
plt.xlim(0, 512)
plt.grid()
plt.title('Amplitude Modulation')
plt.show()
