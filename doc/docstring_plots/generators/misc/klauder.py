#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


import matplotlib.pyplot as plt
from tftb.generators import klauder

x = klauder(128)
plt.plot(x)
plt.xlim(0, 128)
plt.grid()
plt.title('Klauder Wavelet')
plt.show()
