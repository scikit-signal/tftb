#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

import numpy as np
from tftb.generators.api import atoms
from tftb.processing.api import pseudo_wigner_ville
import matplotlib.pyplot as plt

x = np.array([[32, .15, 20, 1],
             [96, .15, 20, 1],
             [32, .35, 20, 1],
             [96, .35, 20, 1]])
g = atoms(128, x)
tfr = pseudo_wigner_ville(g)
threshold = (np.abs(tfr) ** 2) * 0.05
tfr[np.abs(tfr) ** 2 <= threshold] = 0.0

plt.contour(np.abs(tfr) ** 2, levels=range(5))
plt.show()
