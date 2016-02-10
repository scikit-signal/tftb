#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example showing a Bertrand distribution of a hyperbolic group delay signal.
Figure 4.21 from the tutorial.
"""

from tftb.processing.affine import BertrandDistribution
from tftb.generators import gdpower


sig = gdpower(128)[0]
bert = BertrandDistribution(sig, fmin=0.01, fmax=0.22)
bert.run()
bert.plot(kind="contour", show_tf=True)
