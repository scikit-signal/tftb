#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Examples from section 4.1.1 of the tutorial.
"""

import numpy as np
from tftb.generators.api import fmlin
from tftb.processing.api import wigner_ville
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

sig = fmlin(256)[0]
tfr = wigner_ville(sig)
threshold = (np.abs(tfr) ** 2) * 0.05
tfr[np.abs(tfr) ** 2 <= threshold] = 0.0

x = np.arange(256)
y = np.linspace(0, 0.5, 256)
X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, np.abs(tfr), cmap=plt.cm.jet)
ax.set_xlabel('Time')
ax.set_ylabel('Frequency')
ax.set_zlabel('Amplitude')
plt.show()
