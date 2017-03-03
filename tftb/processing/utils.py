#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Miscellaneous processing utilities."""

import numpy as np


def get_spectrum(signal):
    return np.fft.fftshift(np.abs(np.fft.fft(signal)) ** 2)


def integrate_2d(mat, x=None, y=None):
    """integrate_2d

    :param mat:
    :param x:
    :param y:
    :type mat:
    :type x:
    :type y:
:return:
:rtype:
    :Example:
    >>> from __future__ import print_function
    >>> from tftb.generators import altes
    >>> from tftb.processing import Scalogram
    >>> x = altes(256, 0.1, 0.45, 10000)
    >>> tfr, t, f, _ = Scalogram(x).run()
    >>> print("%.3f" % integrate_2d(tfr, t, f))
    2.000
    """
    m, n = mat.shape
    if x is None:
        x = np.arange(n)
    if y is None:
        y = np.arange(m)

    mat = (mat.sum(1) - mat[:, 0] / 2.0 - mat[:, n - 1] / 2.0) * (x[1] - x[0])
    dmat = mat[:(m - 1)] + mat[1:]
    dy = (y[1:] - y[:(m - 1)]) / 2.0
    return np.sum(dmat * dy)


def derive_window(window):
    """Calculate derivative of a window function.

    :param window: Window function to be differentiated. This is expected to be
    a standard window function with an odd length.
    :type window: array-like
    :return: Derivative of the input window
    :rtype: array-like
    :Example:
    >>> from scipy.signal import hanning
    >>> import matplotlib.pyplot as plt                   #doctest: +SKIP
    >>> window = hanning(210)
    >>> derivation = derive_window(window)
    >>> plt.subplot(211), plt.plot(window)                #doctest: +SKIP
    >>> plt.subplot(212), plt.plot(derivation) #doctest: +SKIP

    .. plot:: docstring_plots/processing/utils/derive_window.py
    """
    lh = (window.shape[0] - 1) / 2.0
    step_height = (window[0] + window[-1]) / 2.0
    ramp = (window[0] - window[-1]) / (window.shape[0] - 1)
    base = np.arange(-lh, lh + 1)
    base = window - step_height - ramp * base
    base = np.hstack((np.array([0]), base, np.array([0])))
    dw = (base[2:(window.shape[0] + 2)] - base[:window.shape[0]]) / 2.0 + ramp
    dw[0] += step_height
    dw[-1] -= step_height
    return dw

if __name__ == '__main__':
    x = np.arange(1, 16).reshape(5, 3)
    print(integrate_2d(x))
