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
            ti = timestamps[ti]
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
    from tftb.generators.api import fmlin
    from scipy.signal import kaiser
    from tftb.processing.reassigned import pseudo_wigner_ville
    sig = fmlin(128, 0.1, 0.4)[0]
    fwindow = kaiser(47, beta=3 * np.pi)
    tfr, rtfr, hat = pseudo_wigner_ville(sig, fwindow=fwindow)
    tifd = friedman_density(tfr, hat, 'rpwv')
