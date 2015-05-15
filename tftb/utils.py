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
