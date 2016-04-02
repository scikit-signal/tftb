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
from tftb.processing.utils import integrate_2d


def hough_transform(image, m=None, n=None):
    """hough_transform

    :param image:
    :param m:
    :param n:
    :type image:
    :type m:
    :type n:
    :return:
    :rtype:
    """
    xmax, ymax = image.shape
    if m is None:
        m = xmax
    if n is None:
        n = ymax

    rhomax = np.sqrt((xmax ** 2) + (ymax ** 2)) / 2.0
    deltar = rhomax / (m - 1.0)
    deltat = 2 * np.pi / n

    ht = np.zeros((m, n))
    imax = np.amax(image)

    if xmax % 2 != 0:
        xc = (xmax + 1) // 2
        xf = xc - 1
    else:
        xc = xf = xmax // 2
    x0 = 1 - xc

    if ymax % 2 != 0:
        yc = (ymax + 1) // 2
        yf = yc - 1
    else:
        yc = yf = ymax // 2
    y0 = 1 - yc

    for x in range(x0, xf + 1):
        for y in range(y0, yf + 1):
            if np.abs(image[x + xc - 1, y + yc - 1]) > imax / 20.0:
                for theta in np.linspace(0, 2 * np.pi - deltat, n):
                    rho = x * np.cos(theta) - y * np.sin(theta)
                    if (rho >= 0) and (rho <= rhomax):
                        ht[int(np.round(rho / deltar)),
                           int(np.round(theta / deltat))] += image[x + xc - 1,
                                                                   y + yc - 1]
    rho = np.linspace(0, rhomax, n)
    theta = np.linspace(0, 2 * np.pi - deltat, n)
    return ht, rho, theta


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
        timestamps = np.arange(n) + 1
    elif np.allclose(timestamps, np.arange(n)):
        timestamps += 1
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
    ifrow = iflaws.shape[0]
    if timestamps is None:
        timestamps = np.arange(iflaws[0, :].shape[0])
    if n_fbins is None:
        n_fbins = iflaws[0, :].shape[0]

    tcol = timestamps.shape[0]

    tfr = np.zeros((n_fbins, tcol))
    for icol in range(tcol):
        ti = timestamps[icol]
        for fi in range(ifrow):
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
    if timestamps is None:
        timestamps = np.arange(tfrcol)

    tifd = np.zeros((tfrrow, tfrcol))
    bins = 0.5 + np.arange(tfrrow)
    bin_edges = np.r_[-np.Inf, 0.5 * (bins[:-1] + bins[1:]), np.Inf]
    threshold = np.sum(np.sum(tfr)) * 0.5 / tfr.size

    for j in range(tfrcol):
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
    if timestamps is None:
        timestamps = np.arange(tfrcol)

    n_fbins = tfrrow
    freqs = np.arange(n_fbins)
    threshold = np.sum(np.sum(tfr)) * 0.5 / tfr.size

    if method in ('rpwv', 'rpmh'):
        for icol in range(tfrcol):
            ti = timestamps[icol]
            indices = np.logical_and(tfr[:, icol] > threshold,
                                     re_mat[:, icol] - freqs == 0)
            if np.any(indices):
                time_points = np.ones((indices.sum(),)) * ti
                freq_points = np.arange(indices.shape[0])[indices] / (2.0 * n_fbins)
    elif method == "rspwv":
        for icol in range(tfrcol):
            ti = timestamps[icol]
            condt1 = np.real(re_mat[:, icol]) - freqs == 0
            condt2 = np.imag(re_mat[:, icol]) - icol == 0
            condt3 = tfr[:, icol] > threshold
            indices = np.logical_and(condt1, condt2, condt3)
            if np.any(indices):
                time_points = np.ones((indices.sum(),)) * ti
                freq_points = np.arange(indices.shape[0])[indices] / (2.0 * n_fbins)
    elif method in ("rsp", "type1"):
        for icol in range(tfrcol):
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
    from tftb.generators import atoms
    from tftb.processing import Spectrogram
    s = atoms(64, np.array([[32, 0.3, 16, 1]]))
    spec = Spectrogram(s)
    tfr, t, f = spec.run()
    print(renyi_information(tfr, t, f, 3))
