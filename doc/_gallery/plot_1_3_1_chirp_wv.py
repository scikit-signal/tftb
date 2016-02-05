#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
====================================
Wigner-Ville Distribution of a Chirp
====================================

Construct a chirp signal and visualize its `Wigner-Ville distribution
<https://en.wikipedia.org/wiki/Wigner_distribution_function>`_.

Figure 1.3 from the tutorial.
"""

from tftb.generators import fmlin
from tftb.processing.cohen import WignerVilleDistribution

n_points = 128
fmin, fmax = 0.0, 0.5
signal, _ = fmlin(n_points, fmin, fmax)

# Wigner-Ville distribution of the chirp.

wvd = WignerVilleDistribution(signal)
wvd.run()
wvd.plot(kind='contour', extent=[0, n_points, fmin, fmax])
