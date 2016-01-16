#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
=====================================================
Wigner-Ville Distribution of a Dirac Impulse Function
=====================================================

This example demonstrates the Wigner-Ville distribution of a Dirac impulse
function, and shows the limitations of the WV distribution when applied to
broadband signals.

"""

from tftb.generators import anapulse
from tftb.processing import WignerVilleDistribution
import matplotlib.pyplot as plt
import numpy as np

sig = anapulse(128)
tfr = WignerVilleDistribution(sig).run()[0]
t = np.arange(tfr.shape[1])
f = np.linspace(0, 0.5, tfr.shape[0])

threshold = (np.abs(tfr) ** 2) * 0.05
tfr[(np.abs(tfr) ** 2) <= threshold] = 0.0

plt.figure()
plt.contour(t, f, np.abs(tfr) ** 2, 4)
plt.grid(True)
plt.title('Wigner Ville Distribution of a Dirac impulse.')
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')
plt.show()
