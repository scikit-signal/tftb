#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for the generators.noise module."""

import unittest
import numpy as np
from tftb.tests.test_base import TestBase
from tftb.generators import noise


class TestNoise(TestBase):

    def test_noisecu(self):
        """Test uniform white noise generation."""
        x = noise.noisecu(128)
        self.assertAlmostEqual(x.std() ** 2, 1, places=1)

    def test_noisecg(self):
        """Test Gaussian white noise generation."""
        x = noise.noisecg(128)
        self.assertAlmostEqual(x.std() ** 2, 1, places=1)

    def test_dopnoise(self):
        """Test doppler noise generation"""
        signal, iflaw = noise.dopnoise(500, 200, 60, 10, 70, 128)
        energy = np.sum(np.abs(signal) ** 2)
        self.assertAlmostEqual(energy, 1, 3)
        self.assert_is_monotonic_decreasing(iflaw)
        signal, iflaw = noise.dopnoise(500, 200, 60, 10, 70)
        energy = np.sum(np.abs(signal) ** 2)
        self.assertAlmostEqual(energy, 1, 3)
        self.assert_is_monotonic_decreasing(iflaw)


if __name__ == '__main__':
    unittest.main()
