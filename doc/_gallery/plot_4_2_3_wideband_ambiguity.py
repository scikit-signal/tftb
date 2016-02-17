#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Cube26 product code
#
# (C) Copyright 2015 Cube26 Software Pvt Ltd
# All right reserved.
#
# This file is confidential and NOT open source.  Do not distribute.
#

"""
==============================================
Wideband Ambiguity Function of an Altes Signal
==============================================

For wideband signals, the narrow band ambiguity function is not appropriate for
wideband signals. So we consider a wide band ambiguity function. This function
is equivalent to the wavelet transform of a signal whose mother wavelet is the
signal itself.

Figure 4.25 from the tutorial.

"""

from tftb.generators import altes
from tftb.processing.ambiguity import wide_band
import matplotlib.pyplot as plt
import numpy as np

signal = altes(128, 0.1, 0.45)
waf, tau, theta = wide_band(signal, 0.1, 0.35, 64)
plt.contour(tau, theta, np.abs(waf) ** 2)
plt.xlabel("Delay")
plt.ylabel("Log(Scale)")
plt.title("Wide Band Ambiguity Function")
plt.show()
