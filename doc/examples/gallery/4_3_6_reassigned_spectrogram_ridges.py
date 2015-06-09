#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example of Friedman's instantaneous frequency density calculation.
"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators.api import fmlin, fmsin, fmconst
from tftb.processing.reassigned import spectrogram
from tftb.processing.postprocessing import ridges

sig1, if1 = fmsin(60, 0.16, 0.35, 50, 1, 0.35, 1)
sig2, if2 = fmlin(60, 0.3, 0.1)
sig3, if3 = fmconst(60, 0.4)
sig = np.hstack((sig1, np.zeros((8,)), sig2 + sig3))

tfr, rtfr, hat = spectrogram(sig)
x, y = ridges(tfr, hat)

plt.plot(x, y, '.')
plt.title("Ridges in Spectrogram")
plt.xlabel('Time')
plt.ylabel('Frequency')
plt.show()
