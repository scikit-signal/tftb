#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Wigner Ville distribution of two simultaneous chirps.
"""

from tftb.generators import fmlin, sigmerge
from tftb.processing.cohen import WignerVilleDistribution

N = 64
sig = sigmerge(fmlin(N, 0, 0.4)[0], fmlin(N, 0.3, 0.5)[0], 1)
tfr = WignerVilleDistribution(sig)
tfr.run()
tfr.plot(kind='contour', sqmod=True, show_tf=True)
