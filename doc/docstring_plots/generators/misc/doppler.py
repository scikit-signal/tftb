#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import doppler
import numpy as np
import matplotlib.pyplot as plt

fm, am, iflaw = doppler(512, 200.0, 65.0, 10.0, 50.0)
plt.subplot(211), plt.plot(np.real(am * fm))
plt.title('Doppler')
plt.grid()
plt.xlim(0, 512)
plt.subplot(212), plt.plot(iflaw)
plt.title('Instantaneous Freqeuncy')
plt.grid()
plt.xlim(0, 512)

plt.show()
