#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import amtriang
import matplotlib.pyplot as plt

x = amtriang(160)
plt.plot(x)
plt.title('Triangular Amplitude Modulation')
plt.grid()
plt.show()
