#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import noisecu
import numpy as np
import matplotlib.pyplot as plt

noise = noisecu(512)
print(np.real(noise.mean()))
print(noise.std() ** 2)

plt.subplot(211), plt.plot(np.real(noise))
plt.title('Analytic complex white noise')
plt.xlim(0, 512)
plt.grid()
plt.subplot(212), plt.plot(np.linspace(-0.5, 0.5, 512),
                           np.abs(np.fft.fftshift(np.fft.fft(noise))) ** 2)
plt.title('Energy spectrum')
plt.xlabel('Normalized Frequency')
plt.xlim(-0.5, 0.5)
plt.grid()
plt.show()
