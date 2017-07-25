#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
==========================================================================
Pseudo-Wigner-Ville Distribution of a Gaussian Atom and a Complex Sinusoid
==========================================================================

This example demonstrates the pseudo Wigner Ville distribution of a signal
composed from a Gaussian atom and a complex sinusoid with constant frequency
modulation. Note that the frequency resolution is relatively worse than that of
the Wigner-Ville representation, and the interferences have not been resolved
properly.

Figure 4.9 from the tutorial.
"""

from tftb.generators import fmconst, amgauss
from tftb.processing import PseudoWignerVilleDistribution
import numpy as np

t = np.linspace(0, 1, 128)
sig = fmconst(128, 0.15)[0] + amgauss(128) * fmconst(128, 0.4)[0]
tfr = PseudoWignerVilleDistribution(sig, timestamps=t)
tfr.run()
tfr.plot(show_tf=True, kind="contour")
