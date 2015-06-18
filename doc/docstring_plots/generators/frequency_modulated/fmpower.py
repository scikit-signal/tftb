#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import fmpower
import numpy as np
import matplotlib.pyplot as plt

z, iflaw = fmpower(128, 0.5, (1, 0.5, 100, 0.1))
plt.subplot(211), plt.plot(np.real(z))
plt.xlim(0, 128)
plt.grid()
plt.title('Power Law Modulation')
plt.subplot(212), plt.plot(iflaw)
plt.xlim(0, 128)
plt.grid()
plt.title('Instantaneous Frequency')
plt.show()
