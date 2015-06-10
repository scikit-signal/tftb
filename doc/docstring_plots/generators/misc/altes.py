#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators.misc import altes
import matplotlib.pyplot as plt

x = altes(128, 0.1, 0.45)
plt.plot(x)
plt.xlim(0, 128)
plt.grid()
plt.title("Altes signal")
plt.show()
