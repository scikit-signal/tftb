#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

from tftb.generators import fmlin, sigmerge, noisecg
from tftb.processing.cohen import wigner_ville
import matplotlib.pyplot as plt
import numpy as np

# Generate a chirp signal

n_points = 128
fmin, fmax = 0.0, 0.5

signal, _ = fmlin(n_points, fmin, fmax)

# Noisy chirp

noisy_signal = sigmerge(signal, noisecg(128), 0)


# Wigner-Ville spectrum of noisy chirp.

tfr = wigner_ville(noisy_signal)
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0
plt.contour(tfr, extent=[0, n_points, fmin, fmax])
plt.grid()
plt.title('Wigner-Ville distribution of Noisy Signal (threshold = 5%)')
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')
plt.show()
