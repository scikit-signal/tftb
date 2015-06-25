#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 4.1.4 of the tutorial.
"""

import numpy as np
from tftb.generators import atoms
from tftb.processing import MargenauHillDistribution


sig = atoms(128, np.array([[32, 0.15, 20, 1], [96, 0.32, 20, 1]]))
tfr = MargenauHillDistribution(sig)
tfr.run()
tfr.plot(show_tf=True, kind='contour', sqmod=False, threshold=0,
         contour_y=np.linspace(0, 0.5, tfr.tfr.shape[0] / 2))
