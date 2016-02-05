#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
============================================================
Margenau-Hill Representation of Chirps with Different Slopes
============================================================

This example demonstrates the Margenau-Hill distribution of a signal
composed of two chirps with Gaussian amplitude modulation but havind linear
frequency modulations with different slopes. This distribution too, like the
Wigner-Ville distribution, spearates the signal terms, but produces
interferences such that they appear diagonally opposite their corresponding
signals.

Figure 4.14 from the tutorial.
"""

import numpy as np
from tftb.generators import atoms
from tftb.processing import MargenauHillDistribution


sig = atoms(128, np.array([[32, 0.15, 20, 1], [96, 0.32, 20, 1]]))
tfr = MargenauHillDistribution(sig)
tfr.run()
tfr.plot(show_tf=True, kind='contour', sqmod=False, threshold=0,
         contour_y=np.linspace(0, 0.5, tfr.tfr.shape[0] / 2))
