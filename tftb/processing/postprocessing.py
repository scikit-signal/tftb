#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Postprocessing functions.
"""

import numpy as np
from tftb.utils import init_default_args
from tftb.processing.utils import integrate_2d


def renyi_information(tfr, timestamps=None, freq=None, alpha=3.0):
    """renyi_information

    :param tfr:
    :param timestamps:
    :param freq:
    :param alpha:
    :type tfr:
    :type timestamps:
    :type freq:
    :type alpha:
:return:
:rtype:
    """
    if alpha == 1 and tfr.min().min() < 0:
        raise ValueError("Distribution with negative values not allowed.")
    m, n = tfr.shape
    if timestamps is None:
        timestamps = np.arange(n)
    if freq is None:
        freq = np.arange(m)
    freq.sort()
    tfr = tfr / integrate_2d(tfr, timestamps, freq)
    if alpha == 1:
        R = -integrate_2d(tfr * np.log2(tfr + np.spacing(1)), timestamps, freq)
    else:
        R = np.log2(integrate_2d(tfr ** alpha, timestamps, freq) + np.spacing(1))
        R = R / (1 - alpha)
    return R


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
    timestamps, n_fbins = init_default_args(iflaws[0, :],
            timestamps=timestamps, n_fbins=n_fbins)

    tcol = timestamps.shape[0]

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


def friedman_density(tfr, re_mat, timestamps=None):
    """friedman_density

    :param tfr:
    :param re_mat:
    :param timestamps:
    :type tfr:
    :type re_mat:
    :type timestamps:
:return:
:rtype:
    """
    tfrrow, tfrcol = tfr.shape
    hatrow, hatcol = re_mat.shape
    if timestamps is None:
        timestamps = np.arange(tfrcol)

    tifd = np.zeros((tfrrow, tfrcol))
    bins = 0.5 + np.arange(tfrrow)
    bin_edges = np.r_[-np.Inf, 0.5 * (bins[:-1] + bins[1:]), np.Inf]
    threshold = np.sum(np.sum(tfr)) * 0.5 / tfr.size

    for j in xrange(tfrcol):
        indices = tfr[:, j] > threshold
        if np.any(indices):
            occurences, _ = np.histogram(np.real(re_mat[indices, j]), bin_edges)
            tifd[:, j] = occurences
    tifd = tifd / np.sum(np.sum(tifd))
    return tifd


def ridges(tfr, re_mat, timestamps=None, method='rsp'):
    """ridges

    :param tfr:
    :param re_mat:
    :param timestamps:
    :param method:
    :type tfr:
    :type re_mat:
    :type timestamps:
    :type method:
:return:
:rtype:
    """
    method = method.lower()
    tfrrow, tfrcol = tfr.shape
    hatrow, hatcol = re_mat.shape
    if timestamps is None:
        timestamps = np.arange(tfrcol)

    n_fbins = tfrrow
    freqs = np.arange(n_fbins)
    threshold = np.sum(np.sum(tfr)) * 0.5 / tfr.size

    if method in ('rpwv', 'rpmh'):
        for icol in xrange(tfrcol):
            ti = timestamps[icol]
            indices = np.logical_and(tfr[:, icol] > threshold,
                                     re_mat[:, icol] - freqs == 0)
            if np.any(indices):
                time_points = np.ones((indices.sum(),)) * ti
                freq_points = np.arange(indices.shape[0])[indices] / (2.0 * n_fbins)
    elif method == "rspwv":
        for icol in xrange(tfrcol):
            ti = timestamps[icol]
            condt1 = np.real(re_mat[:, icol]) - freqs == 0
            condt2 = np.imag(re_mat[:, icol]) - icol == 0
            condt3 = tfr[:, icol] > threshold
            indices = np.logical_and(condt1, condt2, condt3)
            if np.any(indices):
                time_points = np.ones((indices.sum(),)) * ti
                freq_points = np.arange(indices.shape[0])[indices] / (2.0 * n_fbins)
    elif method in ("rsp", "type1"):
        for icol in xrange(tfrcol):
            ti = timestamps[icol]
            condt1 = np.real(re_mat[:, icol]) - freqs == 0
            condt2 = np.imag(re_mat[:, icol]) - icol == 0
            condt3 = tfr[:, icol] > threshold
            indices = np.logical_and(condt1, condt2, condt3)
            if np.any(indices):
                time_points = np.ones((indices.sum(),)) * ti
                freq_points = np.arange(indices.shape[0])[indices] / n_fbins
    else:
        raise ValueError("Unknown time frequency representation.")
    return time_points, freq_points

if __name__ == '__main__':
    from tftb.generators.api import atoms
    from tftb.processing.cohen import spectrogram
    s = atoms(64, np.array([[16, .2, 10, 1], [40, 0.4, 12, 1]]))
    tfr, t, f = spectrogram(s)
    print renyi_information(tfr, t, f)
