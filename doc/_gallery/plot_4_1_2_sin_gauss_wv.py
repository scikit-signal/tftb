#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
===================================================================
Wigner Ville distribution of a Gaussian Atom and a Complex Sinusoid
===================================================================

This example demonstrates the Wigner Ville distribution of a signal
composed from a Gaussian atom and a complex sinusoid with constant frequency
modulation. Although the representation does isolate the atom and the sinusoid
as independent phenomena in the signal, it also produces some interference
between them.

Figure 4.8 from the tutorial.
"""

from tftb.generators import fmconst, amgauss
from tftb.processing import WignerVilleDistribution

sig = fmconst(128, 0.15)[0] + amgauss(128) * fmconst(128, 0.4)[0]
tfr = WignerVilleDistribution(sig)
tfr.run()
tfr.plot(show_tf=True, kind='contour')
