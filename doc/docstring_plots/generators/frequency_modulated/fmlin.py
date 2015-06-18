#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import amgauss, fmlin
import numpy as np
import matplotlib.pyplot as plt

z = amgauss(128, 50.0, 30.0) * fmlin(128, 0.05, 0.3, 50)[0]
plt.plot(np.real(z))
plt.xlim(0, 128)
plt.grid()
plt.title('Linear Frequency Modulation')
plt.show()
