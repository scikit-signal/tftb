#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import anapulse
import matplotlib.pyplot as plt
import numpy as np

x = 2.5 * anapulse(512, 301)
plt.plot(np.real(x))
plt.xlim(0, 512)
plt.grid()
plt.title('Analytic Dirac Impulse')
plt.show()
