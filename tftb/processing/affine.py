#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Bilinear Time-Frequency Processing in the Affine Class.
"""

import numpy as np


def bertrand(signal, timestamps=None, fmin=None, fmax=None, n_voices=None):
    """bertrand

    :param signal:
    :param timestamps:
    :param fmin:
    :param fmax:
    :param n_voices:
    :type signal:
    :type timestamps:
    :type fmin:
    :type fmax:
    :type n_voices:
:return:
:rtype:
    """
    xrow = signal.shape[0]
    if timestamps is None:
        timestamps = np.arange(xrow)

    tcol = timestamps.shape[0]
    x1 = signal.copy()
    x2 = signal.copy()

    s1 = np.real(x1)
    s2 = np.real(x2)
    m = (xrow + (xrow % 2)) / 2
    t = np.arange(xrow) - m - 1
    tmin = 1
    tmax = xrow
    T = tmax - tmin

    if fmin is None:
        stf1 = np.fft.fft(np.fft.fftshift(s1[timestamps.min():(timestamps.max() + 1)]))
        nstf = stf1.shape[0]
        sp1 = np.abs(stf1[:int(np.floor(nstf / 2.0))]) ** 2
        f = np.linspace(0, 0.05, np.floor(nstf / 2.0) + 1)[:int(np.floor(nstf / 2.0))]
        fmin = max([0.01, 0.05 * np.floor(f[indmin] / 0.05)])

