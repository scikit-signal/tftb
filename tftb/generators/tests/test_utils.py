#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Cube26 product code
#
# (C) Copyright 2015 Cube26 Software Pvt Ltd
# All right reserved.
#
# This file is confidential and NOT open source.  Do not distribute.
#

"""
Tests for tftb.generators.utils
"""

import unittest
import numpy as np

from tftb.generators import utils, fmlin


class TestUtils(unittest.TestCase):

    def test_sigmerge(self):
        """Test merging of signals with a given SNR."""
        signal = fmlin(128)[0]
        noise = np.random.randn(128,)
        noisy_signal = utils.sigmerge(signal, noise)
        gamma_estimate = np.sqrt(signal.var() / noise.var())
        np.testing.assert_allclose(noisy_signal,
                                   signal + gamma_estimate * noise, rtol=1e-2,
                                   atol=1e-2)


if __name__ == '__main__':
    unittest.main()
