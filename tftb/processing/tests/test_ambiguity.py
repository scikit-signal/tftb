#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Tests for tftb.processing.ambiguity
"""

import unittest
import numpy as np

from tftb.generators import fmlin, amgauss, altes
from tftb.processing import ambiguity
from tftb.tests.test_base import TestBase


class TestWideBand(TestBase):
    """Tests for the wide band ambiguity function."""
    def setUp(self):
        self.signal = altes(128, 0.1, 0.45)
        self.tfr, self.tau, self.theta = ambiguity.wide_band(self.signal)

    def test_max_value(self):
        """Test if the maximum value of the ambiguity function occurs at the
        origin"""
        max_ix = np.argmax(np.abs(self.tfr))
        self.assertAlmostEqual(float(self.tfr.size) / max_ix, 2, places=1)


class TestNarrowBand(TestBase):
    """Tests for the narrow band ambiguity function"""

    def setUp(self):
        x = fmlin(64, 0.2, 0.5)[0] * amgauss(64)
        y = fmlin(64, 0.3, 0)[0] * amgauss(64)
        self.signal = np.hstack((x, y))
        self.tfr, self.lag, self.doppler = ambiguity.narrow_band(self.signal)

    def test_max_value(self):
        """Test if the maximum value of the ambiguity function occurs at the
        origin"""
        xorg, yorg = map(lambda x: int(x / 2), self.tfr.shape)
        abs_mat = np.abs(self.tfr) ** 2
        max_val = abs_mat[xorg, yorg]
        for i in range(abs_mat.shape[0]):
            for j in range(abs_mat.shape[1]):
                self.assertTrue(abs_mat[i, j] <= max_val)

    def test_volume_invariance(self):
        """Test the volume invariance property of the narrow band ambiguity
        function."""
        volume = np.abs(self.tfr[self.doppler == 0, self.lag == 0]) ** 2
        volume_integral = (np.abs(self.tfr) ** 2).sum().sum()
        self.assertAlmostEqual(volume[0],
                               volume_integral / (self.signal.shape[0] / 2))

    @unittest.skip("Not quite ready yet.")
    def test_symmetry(self):
        """Test the symmetry property of the narrow band ambiguity function."""
        tfr = self.tfr[1:, :]
        self.assert_is_hermitian(tfr)


if __name__ == '__main__':
    unittest.main()
