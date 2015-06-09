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
from tftb.generators.api import atoms
from tftb.processing.api import renyi_information
from tftb.processing.cohen import wigner_ville

# FIXME: This is fully wrong.

sig = atoms(128, np.array([[64, 0.25, 20, 1]]))
tfr = wigner_ville(sig)
print renyi_information(tfr)  # -0.2075

sig = atoms(128, np.array([[32, 0.25, 20, 1], [96, 0.25, 20, 1]]))
tfr = wigner_ville(sig)

print renyi_information(tfr)  # 0.77

sig = atoms(128, np.array([[32, 0.15, 20, 1], [96, 0.15, 20, 1],
                           [32, 0.35, 20, 1], [96, 0.35, 20, 1]]))
tfr = wigner_ville(sig)
print renyi_information(tfr)  # 1.8029
