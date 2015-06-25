#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Examples from section 3.4.1 of the tutorial.
"""

from tftb.generators import fmlin
from tftb.processing.cohen import Spectrogram
import numpy as np
import matplotlib.pyplot as plt

sig = fmlin(128, 0, 0.3)[0] + fmlin(128, 0.2, 0.5)[0]
window = np.exp(np.log(0.005) * np.linspace(-1, 1, 23) ** 2)
spec = Spectrogram(sig, fwindow=window, n_fbins=128)
spec.run()
spec.plot(show_tf=True, cmap=plt.cm.gray)
