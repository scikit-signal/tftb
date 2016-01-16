#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from seciton 4.1.1 of the tutorial.
"""

import numpy as np
from tftb.generators import atoms
from tftb.processing import WignerVilleDistribution
import matplotlib.pyplot as plt

x = np.array([[32, .15, 20, 1],
             [96, .32, 20, 1]])
g = atoms(128, x)
tfr = WignerVilleDistribution(np.real(g)).run()[0]
threshold = (np.abs(tfr) ** 2) * 0.05
tfr[np.abs(tfr) ** 2 <= threshold] = 0.0

plt.contour(np.abs(tfr) ** 2, levels=list(range(5)))
plt.show()
