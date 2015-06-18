#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators import dopnoise
from tftb.processing.freq_domain import inst_freq

z, iflaw = dopnoise(500, 200.0, 60.0, 10.0, 70.0, 128.0)
plt.subplot(211), plt.plot(np.real(z))
plt.grid()
plt.title('Complex noisy Doppler signal')
plt.xlim(0, 500)
ifl, t = inst_freq(z, np.arange(11, 479), 10)
plt.subplot(212)
plt.plot(iflaw, 'r', label='actual')
plt.plot(t, ifl, 'g', label='estimated')
plt.xlim(0, 500)
plt.legend()
plt.grid()
plt.title('Instantaneous frequency')
plt.show()
