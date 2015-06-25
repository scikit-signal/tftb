#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 1.3.1
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
