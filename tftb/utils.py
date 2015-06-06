#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Miscellaneous utilities."""

import numpy as np
import warnings


def init_default_args(signal, **kwargs):
    """init_default_args

    :param x:
    :param **kwargs:
    :type x:
    :type **kwargs:
:return:
:rtype:
    """
    for varname, value in kwargs.iteritems():
        if "time" in varname:
            if value is None:
                kwargs[varname] = np.arange(signal.shape[0])
        elif ("n_fbins" in varname) or ("freq_bins" in varname):
            if value is None:
                kwargs[varname] = signal.shape[0]
            elif 2 ** nextpow2(value) != value:
                msg = "For faster computations, n_fbins should be a power of 2."
                warnings.warn(msg, UserWarning)

    return kwargs.values()


def izak(x):
    """Inverse Zak transform."""
    if x.ndim == 2:
        n, m = x.shape
    else:
        n, m = x.shape[0], 1
    sig = np.zeros((n * m, ), dtype=complex)
    for im in xrange(m):
        sig[im + np.arange(n) * m] = np.sqrt(n) * np.fft.ifft(x[:, im], axis=0)
    return sig


def nextpow2(n):
    """Compute the exponent of the next higher power of 2.

    :param n: Number whose next higest power of 2 needs to be computed.
    :type n: int, np.ndarray
    :rtype: int, np.ndarray
    """
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return m_i


def divider(N):
    """Compute two factors of N such that they are as close as possible to sqrt(N).

    :param N: Number to be divided.
    """
    n = np.floor(np.sqrt(N))
    while True:
        old = n
        m = np.ceil(N / float(n))
        n = np.floor(N / float(m))
        if n == old:
            break
    return n, m


def nearest_odd(N):
    """Get the nearest odd number for each value of N."""
    if isinstance(N, np.ndarray):
        y = np.floor(N)
        y[np.remainder(y, 2) == 0] = np.ceil(N[np.remainder(y, 2) == 0])
        y[np.remainder(y, 2) == 0] += 1
        return y
    if N % 2 == 0:
        return N + 1
    return N


def modulo(x, N):
    """Compute the congruence of each element of x modulo N.

    :type x: array-like
    :type N: int
    :return: array-like
    """
    if any(np.isreal(x)):
        y = np.mod(x, N)
        y[y == 0] = N
    else:
        y = np.mod(np.real(x), N) + 1j * np.mod(np.imag(x), N)
    return y


if __name__ == '__main__':
    print init_default_args(np.arange(10), timestamps=range(10), n_fbins=3)
