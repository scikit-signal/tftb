#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators import amgauss, fmlin
from tftb.processing import group_delay

x = amgauss(128, 64.0, 30) * fmlin(128, 0.1, 0.4)[0]
fnorm = np.arange(0.1, 0.38, step=0.04)
gd = group_delay(x, fnorm)
plt.plot(gd, fnorm)
plt.xlim(0, 128)
plt.grid()
plt.title('Group delay estimation of linear chirp')
plt.xlabel('Group delay')
plt.ylabel('Normalized frequency')
plt.show()
