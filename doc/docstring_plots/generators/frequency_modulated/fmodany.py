#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.
"""

"""

from tftb.generators import fmlin, fmodany, fmsin
import numpy as np
import matplotlib.pyplot as plt

y1, ifl1 = fmlin(100)
y2, ifl2 = fmsin(100)
iflaw = np.append(ifl1, ifl2)
sig = fmodany(iflaw)

plt.subplot(211), plt.plot(np.real(sig))
plt.grid()
plt.title('Linear and Sinusoidal modulated signal')
plt.subplot(212), plt.plot(iflaw)
plt.grid()
plt.title('Instantaneous frequency')
plt.show()
