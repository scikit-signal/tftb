#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Miscellaneous utilities."""

import numpy as np


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


if __name__ == '__main__':
    print nearest_odd(np.arange(10))
