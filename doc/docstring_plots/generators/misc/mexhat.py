#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import mexhat
import matplotlib.pyplot as plt

plt.plot(mexhat())
plt.grid()
plt.title('Mexican Hat Wavelet')
plt.show()
