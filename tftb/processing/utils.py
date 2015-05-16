#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Miscellaneous processing utilities."""

import numpy as np


def derive_window(window):
    """Calculate derivative of a window function.

    :param window: Window function to be differentiated. This is expected to be
    a standard window function with an odd length.
    :type window: array-like
    :return: Derivative of the input window
    :rtype: array-like
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
    import matplotlib.pyplot as plt
    from scipy.signal import hanning
    window = hanning(210)
    derivative = derive_window(window)
    plt.plot(window, 'b')
    plt.plot(derivative, 'g')
    plt.show()
