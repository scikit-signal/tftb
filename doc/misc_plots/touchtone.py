#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Cube26 product code
#
# (C) Copyright 2015 Cube26 Software Pvt Ltd
# All right reserved.
#
# This file is confidential and NOT open source.  Do not distribute.
#

"""

"""
import numpy as np
import matplotlib.pyplot as plt
fs = 32768
ts = np.linspace(0, 1, fs)
y1 = np.sin(2 * np.pi * 697 * ts)
y2 = np.sin(2 * np.pi * 1336 * ts)
y = (y1 + y2) / 2
plt.plot(ts, y)
plt.xlim(0, 0.1)
plt.show()
