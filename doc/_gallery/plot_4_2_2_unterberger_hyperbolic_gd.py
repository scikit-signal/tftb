#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
===========================================================
Unterberger distribution of a hyperbolic group delay signal
===========================================================

The active Unterberger distribution is the only localized bi-frequency kernel
distribution which localizes perfectly signals having a group delay in
:math:`1/\\nu^{2}`

Figure 4.23 from the tutorial.
"""

from tftb.processing import UnterbergerDistribution
from tftb.generators import gdpower

sig = gdpower(128, -1)[0]
dist = UnterbergerDistribution(sig, fmin=0.01, fmax=0.22, n_voices=172)
dist.run()
dist.plot()
