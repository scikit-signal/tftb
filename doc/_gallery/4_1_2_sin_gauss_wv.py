#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 4.1.2 of the tutorials.
"""

from tftb.generators import fmconst, amgauss
from tftb.processing import WignerVilleDistribution
import numpy as np

sig = fmconst(128, 0.15)[0] + amgauss(128) * fmconst(128, 0.4)[0]
tfr = WignerVilleDistribution(sig)
tfr.run()
tfr.plot(show_tf=True, kind='contour',
        freq_x=(abs(np.fft.fftshift(np.fft.fft(sig))) ** 2)[::-1][:64],
        freq_y=np.arange(sig.shape[0] / 2))
