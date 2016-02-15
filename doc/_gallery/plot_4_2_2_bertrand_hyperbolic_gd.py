#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
========================================================
Bertrand Distribution of a Hyperbolic Group Delay Signal
========================================================

This example shows the Bertrand distribution of a signal with hyperbolic group
delay. The distribution is well localized around the hyperbola, but not
perfectly. The Bertrand distribution operates only on a part of the frequency
range between two bounds :math:`f_{min}` and :math:`f_{max}`.

Figure 4.21 from the tutorial.
"""

from tftb.processing.affine import BertrandDistribution
from tftb.generators import gdpower


sig = gdpower(128)[0]
bert = BertrandDistribution(sig, fmin=0.01, fmax=0.22)
bert.run()
bert.plot()
