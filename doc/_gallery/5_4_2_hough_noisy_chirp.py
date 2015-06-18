#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example showing use of Hough transform on a Wigner-Ville distribution.
"""

import numpy as np
from tftb.generators import noisecg, sigmerge, fmlin
from tftb.processing.cohen import wigner_ville
from tftb.processing.postprocessing import hough_transform
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

N = 64
sig = sigmerge(fmlin(N, 0, 0.3)[0], noisecg(N), 1)
tfr = wigner_ville(sig)

ht, rho, theta = hough_transform(tfr, N, N)
theta, rho = np.meshgrid(theta, rho)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_wireframe(theta, rho, ht)
ax.set_xlabel('Theta')
ax.set_ylabel('Rho')
plt.show()
