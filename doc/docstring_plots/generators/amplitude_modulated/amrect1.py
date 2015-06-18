#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import amrect
import matplotlib.pyplot as plt

x = amrect(160, 90, 40.0)
plt.plot(x)
plt.grid()
plt.title('Rectangular Amplitude Modulation.')
plt.show()

