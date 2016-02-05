#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
=================================
Group Delay Estimation of a Chirp
=================================

Constuct a chirp and estimates its `group delay
<https://en.wikipedia.org/wiki/Group_delay_and_phase_delay>`_.

Figure 2.4 from the tutorial.
"""

from tftb.generators import fmlin
from tftb.processing import group_delay
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
