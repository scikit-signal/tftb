#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
================================================================================
Sampling Effects on the Wigner-Ville Distribution of a Real Valued Gaussian Atom
================================================================================

This example shows the Wigner-Ville distribution of a real valued Gaussian
atom. If a signal is sampled at the Nyquist rate, the WVD is affected by
spectral aliasing and many additional interferences. To fix this, either the
signal may be oversampled, or an analytical signal may be used.

Figure 4.6 from the tutorial.
"""

import numpy as np
from tftb.generators import atoms
from tftb.processing import WignerVilleDistribution

x = np.array([[32, .15, 20, 1],
             [96, .32, 20, 1]])
g = atoms(128, x)
spec = WignerVilleDistribution(np.real(g))
spec.run()
spec.plot(kind="contour", show_tf=True, scale="log")
