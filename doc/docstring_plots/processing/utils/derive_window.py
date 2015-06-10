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
from scipy.signal import hanning
from tftb.processing.utils import derive_window

window = hanning(210)
plt.subplot(211), plt.plot(window)
plt.xlim(0, 210)
plt.grid()
plt.title('Hanning window')
plt.subplot(212), plt.plot(derive_window(window))
plt.xlim(0, 210)
plt.grid()
plt.title('Approximate derivative')
plt.show()
