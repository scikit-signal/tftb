#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Miscellaneous functions.
"""

import numpy as np


def ideal_tfr(iflaws, timestamps=None, n_fbins=None):
    """ideal_tfr

    :param iflaws:
    :param timestamps:
    :param n_fbins:
    :type iflaws:
    :type timestamps:
    :type n_fbins:
:return:
:rtype:
    """
    ifrow, ifcol = iflaws.shape

    if timestamps is None:
        timestamps = np.arange(ifcol)
    tcol = timestamps.shape[0]

    if n_fbins is None:
        n_fbins = ifcol

    tfr = np.zeros((n_fbins, tcol))
    for icol in xrange(tcol):
        ti = timestamps[icol]
        for fi in xrange(ifrow):
            if np.isnan(iflaws[fi, ti]):
                tfr[ti, fi] = np.nan
            else:
                tfr[int(np.round(iflaws[fi, ti] * 2 * (n_fbins - 1))), icol] = 1
    freqs = np.arange(n_fbins, dtype=float) / n_fbins * 0.5
    return tfr, timestamps, freqs

if __name__ == '__main__':
    from tftb.generators.api import fmlin, fmsin
    import matplotlib.pyplot as plt
    n_points = 140
    t = np.arange(n_points)
    x1, if1 = fmlin(n_points, 0.05, 0.3)
    x2, if2 = fmsin(70, 0.35, 0.45, 60)
    if2 = np.hstack((np.zeros((35,)) * np.nan, if2, np.zeros((35,)) * np.nan))
    tfr, t, f = ideal_tfr(np.vstack((if1, if2)))
    plt.contour(t, f, tfr, 1)
    plt.grid(True)
    plt.show()
