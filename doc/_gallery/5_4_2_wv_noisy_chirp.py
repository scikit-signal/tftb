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

from tftb.generators import noisecg, sigmerge, fmlin
from tftb.processing.cohen import WignerVilleDistribution

N = 64
sig = sigmerge(fmlin(N, 0, 0.3)[0], noisecg(N), 1)
wvd = WignerVilleDistribution(sig)
wvd.run()
wvd.plot(kind='contour', show_tf=True, sqmod=True)
