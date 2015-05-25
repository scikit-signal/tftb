#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 1.3.1 of the tutorial.
"""

from tftb.generators.api import fmlin
from tftb.processing.api import group_delay
import numpy as np
import matplotlib.pyplot as plt

signal, _ = fmlin(256)
fnorm = np.linspace(0, .5, 10)
gd = group_delay(signal, fnorm)
plt.plot(gd, fnorm)
plt.grid(True)
plt.xlim(0, 256)
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')
plt.title('Group Delay Estimation')
plt.show()
