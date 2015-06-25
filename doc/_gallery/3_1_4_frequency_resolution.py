#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example from section 3.1.4 of the tutorial.
"""

import numpy as np
from tftb.processing.linear import ShortTimeFourierTransform
from tftb.generators import fmlin, amgauss
import matplotlib.pyplot as plt

x = np.real(amgauss(128) * fmlin(128)[0])
window = np.ones((128,))
stft = ShortTimeFourierTransform(x, n_fbins=128, fwindow=window)
stft.run()
stft.plot(show_tf=True, cmap=plt.cm.gray)
