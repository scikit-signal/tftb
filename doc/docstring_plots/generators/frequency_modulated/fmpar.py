#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import fmpar
import numpy as np
import matplotlib.pyplot as plt

z, iflaw = fmpar(128, (0.4, -0.0112, 8.6806e-05))
plt.subplot(211), plt.plot(np.real(z))
plt.xlim(0, 128)
plt.grid()
plt.title('Parabolic Frequency Modulation')
plt.subplot(212), plt.plot(iflaw)
plt.xlim(0, 128)
plt.grid()
plt.title('Instantaneous Frequency')
plt.show()
