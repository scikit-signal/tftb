#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

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
        gamma = 0.1
        x = utils.sigmerge(signal, noise, gamma)
        h_est = np.linalg.norm(signal) / np.linalg.norm(noise) * 10 ** (-gamma / 20)
        x_hat = signal + h_est * noise
        np.testing.assert_allclose(x, x_hat)

if __name__ == '__main__':
    unittest.main()
