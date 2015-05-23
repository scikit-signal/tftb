#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 1.3.2 of the tutorial.
"""

from os.path import abspath, dirname, join
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from tftb.processing.cohen import pseudo_wigner_ville

DATAPATH = join(abspath(dirname(__file__)), "data", "bat.mat")

bat_signal = loadmat(DATAPATH)['bat'].ravel()

t0 = np.linspace(0, 2500 / 2304.0, 2500)
n = bat_signal.shape[0]
plt.plot(t0[:n], bat_signal)
plt.xlabel('Time [ms]')
plt.title('Sonar Signal from a bat')
plt.grid()
plt.show()

# Energy spectrum of the bat

dsp = np.fft.fftshift(np.abs(np.fft.fft(bat_signal)) ** 2)
f0 = np.arange(-1250, 1250) * 230.4 / 2500
plt.plot(f0[:n], dsp)
plt.xlabel('Normalized Frequency')
plt.ylabel('squared modulus')
plt.title('Energy spectrum of a Bat Sonar Signal')
plt.grid()
plt.show()

# Pseudo Wigner Ville spectrum of the bat signal
tfr = pseudo_wigner_ville(bat_signal)
# threshold = np.amax(tfr) * 0.05
# tfr[tfr <= threshold] = 0.0
plt.contour(tfr)
plt.grid()
plt.title('Wigner-Ville distribution of Noisy Signal (threshold = 5%)')
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')
plt.show()
