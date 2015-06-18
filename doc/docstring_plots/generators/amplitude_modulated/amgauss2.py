#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import amgauss
import matplotlib.pyplot as plt

x = amgauss(160, 90)
plt.plot(x)
plt.grid()
plt.title('Gaussian Amplitude Modulation')
plt.show()
