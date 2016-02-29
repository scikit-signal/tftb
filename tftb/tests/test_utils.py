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
Tests for tftb.utils
"""

import unittest
import numpy as np
from tftb import utils


class TestUtils(unittest.TestCase):

    def test_is_linear(self):
        x = np.arange(10)
        self.assertTrue(utils.is_linear(x))
        x = np.sin(x)
        self.assertFalse(utils.is_linear(x))

    def test_nextpow2(self):
        self.assertEqual(utils.nextpow2(2), 1)
        self.assertEqual(utils.nextpow2(17), 5)
        import warnings
        with warnings.catch_warnings(record=True) as catcher:
            utils.nextpow2(-3)
            self.assertEqual(len(catcher), 1)
            self.assertTrue(catcher[-1].category, RuntimeWarning)

    def test_divider(self):
        self.assertItemsEqual(utils.divider(4), (2, 2))
        self.assertItemsEqual(utils.divider(17), (1, 17))
        self.assertItemsEqual(utils.divider(60), (10, 6))
        x = np.arange(1, 101)
        lowers = np.zeros(x.shape)
        uppers = np.zeros(x.shape)
        for i, num in enumerate(x):
            a, b = utils.divider(num)
            lowers[i] = a
            uppers[i] = b
        perfect_squares = np.arange(1, 11) ** 2
        np.testing.assert_allclose(perfect_squares, x[lowers == uppers])

    def test_nearest_odd(self):
        self.assertEqual(utils.nearest_odd(0), 1)
        self.assertEqual(utils.nearest_odd(2), 3)
        self.assertEqual(utils.nearest_odd(-0.00001), -1)

    def test_modulo(self):
        x = np.arange(1, 11)
        np.testing.assert_allclose(utils.modulo(x, 1), np.ones(x.shape))
        np.testing.assert_allclose(utils.modulo(x, 2),
                np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2]))
        np.testing.assert_allclose(utils.modulo(x, 3),
                np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1]))
        np.testing.assert_allclose(utils.modulo(x, 4),
                np.array([1, 2, 3, 4, 1, 2, 3, 4, 1, 2]))
        np.testing.assert_allclose(utils.modulo(x, 5),
                np.array([1, 2, 3, 4, 5, 1, 2, 3, 4, 5]))

if __name__ == '__main__':
    unittest.main()
