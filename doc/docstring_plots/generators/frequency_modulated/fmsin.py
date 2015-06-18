#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import fmsin
import numpy as np
import matplotlib.pyplot as plt

z = fmsin(140, period=100, t0=20.0, fnorm0=0.3, pm1=-1)[0]
plt.plot(np.real(z))
plt.grid()
plt.title('Sinusoidal Frequency Modulation')
plt.show()
