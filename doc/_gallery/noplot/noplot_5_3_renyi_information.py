#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Examples showing Renyi information measurement.
"""

import numpy as np
from scipy.io import loadmat
from tftb.generators import atoms
from tftb.processing import renyi_information
from tftb.processing.cohen import WignerVilleDistribution

sig = atoms(128, np.array([[64, 0.25, 20, 1]]))
tfr, t, f = WignerVilleDistribution(sig).run()
ideal = loadmat("/tmp/foo.mat")
print(renyi_information(tfr, t, f))  # -0.2075

sig = atoms(128, np.array([[32, 0.25, 20, 1], [96, 0.25, 20, 1]]))
tfr, t, f = WignerVilleDistribution(sig).run()
print(renyi_information(tfr, t, f))  # 0.77

sig = atoms(128, np.array([[32, 0.15, 20, 1], [96, 0.15, 20, 1],
                           [32, 0.35, 20, 1], [96, 0.35, 20, 1]]))
tfr, t, f = WignerVilleDistribution(sig).run()
print(renyi_information(tfr, t, f))  # 1.8029
