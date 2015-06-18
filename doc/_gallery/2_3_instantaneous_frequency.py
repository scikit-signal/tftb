#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 2.3 of the tutorial.
"""

from tftb.generators import fmlin
from tftb.processing.api import plotifl, inst_freq
import numpy as np


signal, _ = fmlin(256)
time_samples = np.arange(3, 257)
ifr = inst_freq(signal)[0]
plotifl(time_samples, ifr)
