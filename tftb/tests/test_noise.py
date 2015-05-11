#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for the generators.noise module."""

import unittest
from tftb.tests.base import TestBase
from tftb.generators import noise


class TestNoise(TestBase):

    def test_noisecu(self):
        x = noise.noisecu(128)
        self.assertAlmostEqual(x.std() ** 2, 1, places=1)

    def test_noisecg(self):
        x = noise.noisecg(128)
        self.assertAlmostEqual(x.std() ** 2, 1, places=1)

if __name__ == '__main__':
    unittest.main()
