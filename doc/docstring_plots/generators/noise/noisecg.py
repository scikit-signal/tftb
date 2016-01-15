#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

import matplotlib.pyplot as plt
import numpy as np
from tftb.generators import noisecg

noise = noisecg(512)
print(noise.mean())
print(noise.std() ** 2)
plt.subplot(211), plt.plot(np.real(noise))
plt.xlim(0, 512)
plt.grid()
plt.title('Analytic complex Gaussian noise.')
plt.subplot(212), plt.plot(np.linspace(-0.5, 0.5, 512),
                           abs(np.fft.fftshift(np.fft.fft(noise))) ** 2)
plt.grid()
plt.title('Energy spectrum')
plt.xlim(-0.5, 0.5)
plt.show()
