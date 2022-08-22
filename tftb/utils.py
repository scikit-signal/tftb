#! /usr/bin/env python

"""Miscellaneous utilities."""

import math
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
    return math.ceil(math.log2(n))


def divider(N):
    """
    Compute two factors of N such that they are as close as possible to sqrt(N).

    :param N: Number to be divided.
    :type N: int
    :return: A tuple of two integers such that their product is `N` and they
        are the closest possible to :math:`\sqrt(N)`  # NOQA: W605
    :rtype: tuple(int)
    :Example:
    >>> from __future__ import print_function
    >>> print(divider(256))
    (16, 16)
    >>> print(divider(10))
    (2, 5)
    >>> print(divider(101))
    (1, 101)
    """
    s = math.sqrt(N)
    if s % 1 == 0:
        s = int(s)
        return s, s
    s = math.ceil(s)
    for i in range(s, 0, -1):
        factor = N // i
        if N % i == 0:
            break
    return i, factor


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
    if isinstance(N, (list, np.ndarray)):
        y = np.floor(N)
        y[np.remainder(y, 2) == 0] += 1
        return y.astype(int)
    if N % 2 == 0:
        return N + 1
    elif math.floor(N) % 2 == 0:
        return math.ceil(N)
    elif math.floor(N) % 2 != 0:
        return math.floor(N)
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
