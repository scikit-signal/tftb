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
    """
    Check if an array is linear.

    :param x: Array to be checked for linearity.
    :param decimals: decimal places upto which the derivative of the array
        should be rounded off (default=5)
    :type x: numpy.ndarray
    :type decimals: int
    :return: If the array is linear
    :rtype: boolean
    :Example:
    >>> import numpy as np
    >>> x = np.linspace(0, 2 * np.pi, 100)
    >>> is_linear(x)
    True
    >>> is_linear(np.sin(x))
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
    :Example:
    >>> from __future__ import print_function
    >>> import numpy as np
    >>> x = np.arange(1, 9)
    >>> print(nextpow2(x))
    [ 0.  1.  2.  2.  3.  3.  3.  3.]
    """
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return m_i


def divider(N):
    """
    Compute two factors of N such that they are as close as possible to sqrt(N).

    :param N: Number to be divided.
    :type N: int
    :return: A tuple of two integers such that their product is `N` and they
        are the closest possible to :math:`\sqrt(N)`
    :rtype: tuple(int)
    :Example:
    >>> from __future__ import print_function
    >>> print(divider(256))
    (16.0, 16.0)
    >>> print(divider(10))
    (2.0, 5.0)
    >>> print(divider(101))
    (1.0, 101.0)
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
    """
    Get the nearest odd number for each value of N.

    :param N: int / sequence of ints
    :return: int / sequence of ints
    :Example:
    >>> from __future__ import print_function
    >>> print(nearest_odd(range(1, 11)))
    [  1.   3.   3.   5.   5.   7.   7.   9.   9.  11.]
    >>> nearest_odd(0)
    1
    >>> nearest_odd(3)
    3.0
    """
    if hasattr(N, "__iter__"):
        N = np.array(N)
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
    """
    Compute the congruence of each element of x modulo N.

    :type x: array-like
    :type N: int
    :return: array-like
    :Example:
    >>> from __future__ import print_function
    >>> print(modulo(range(1, 11), 2))
    [1 2 1 2 1 2 1 2 1 2]
    """
    if any(np.isreal(x)):
        y = np.mod(x, N)
        y[y == 0] = N
    else:
        y = np.mod(np.real(x), N) + 1j * np.mod(np.imag(x), N)
    return y
