#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for the generators.misc module."""

import unittest
from tftb.tests.base import TestBase
from tftb.generators import misc
import numpy as np


class TestMisc(TestBase):

    def test_altest(self):
        ideal = np.array([0.00200822, -0.21928398, 0.66719239, 0.66719239,
                          -0.17666382, -0.17009953, -0.038399, -0.00083597])
        actual = misc.altes(8, 0.1, 0.5)
        np.testing.assert_allclose(ideal, actual, atol=1e-8, rtol=1e-8)

    def test_doppler(self):
        fm, am, iflaw = misc.doppler(512, 200.0, 65, 10, 50)
        self.assert_is_monotonic_decreasing(iflaw)

if __name__ == '__main__':
    unittest.main()
