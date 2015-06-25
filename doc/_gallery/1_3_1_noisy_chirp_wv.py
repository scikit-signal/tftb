#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import fmlin, sigmerge, noisecg
from tftb.processing.cohen import WignerVilleDistribution

# Generate a chirp signal

n_points = 128
fmin, fmax = 0.0, 0.5

signal, _ = fmlin(n_points, fmin, fmax)

# Noisy chirp

noisy_signal = sigmerge(signal, noisecg(128), 0)


# Wigner-Ville spectrum of noisy chirp.

wvd = WignerVilleDistribution(noisy_signal)
wvd.run()
wvd.plot(kind='contour')
