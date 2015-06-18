#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import amexpos
import matplotlib.pyplot as plt

x = amexpos(160)
plt.plot(x)
plt.grid()
plt.title("Two sided exponential amplitude modulation")
plt.show()
