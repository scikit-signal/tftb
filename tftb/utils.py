#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Miscellaneous utilities."""

import numpy as np


def is_linear(x, decimals=5):
    """Check if an array is linear.

    :param x: Array to be checked for linearity.
    :param decimals: decimal places upto which the derivative of the array
        should be rounded off (default=5)
    :type x: numpy.ndarray
    :type decimals: int
    :return: If the array is linear
    :rtype: boolean
    :Example:
    >>> x = np.linspace(0, 2*pi, 100)
    >>> is_linear(x)
    True
    >>> is_linear(sin(x))
    False
    """
    derivative = np.diff(x)
    derivative = np.around(derivative, decimals)
    return np.unique(derivative).shape[0] == 1


def izak(x):
    """Inverse Zak transform."""
    if x.ndim == 2:
        n, m = x.shape
    else:
        n, m = x.shape[0], 1
    sig = np.zeros((n * m, ), dtype=complex)
    for im in range(m):
        sig[im + np.arange(n) * m] = np.sqrt(n) * np.fft.ifft(x[:, im], axis=0)
    return sig


def nextpow2(n):
    """
    Compute the integer exponent of the next higher power of 2.

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
    elif np.floor(N) % 2 == 0:
        return np.ceil(N)
    elif np.floor(N) % 2 != 0:
        return np.floor(N)
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
