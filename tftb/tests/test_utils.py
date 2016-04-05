#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Tests for tftb.utils
"""

import unittest
import six
import numpy as np
from tftb import utils


class TestUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if six.PY3:
            cls.assertItemsEqual = cls.assertSequenceEqual

    def test_is_linear(self):
        """Test the is_linear function."""
        x = np.arange(10)
        self.assertTrue(utils.is_linear(x))
        x = np.sin(x)
        self.assertFalse(utils.is_linear(x))

    def test_nextpow2(self):
        """Test the nextpow2 function."""
        self.assertEqual(utils.nextpow2(2), 1)
        self.assertEqual(utils.nextpow2(17), 5)
        import warnings
        with warnings.catch_warnings(record=True) as catcher:
            utils.nextpow2(-3)
            self.assertEqual(len(catcher), 1)
            self.assertTrue(catcher[-1].category, RuntimeWarning)

    def test_divider(self):
        """Test the divider function."""
        self.assertItemsEqual(utils.divider(4), (2, 2))
        self.assertItemsEqual(utils.divider(17), (1, 17))
        self.assertItemsEqual(utils.divider(60), (6, 10))
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
        """Test the nearest_odd function."""
        self.assertEqual(utils.nearest_odd(0), 1)
        self.assertEqual(utils.nearest_odd(2), 3)
        self.assertEqual(utils.nearest_odd(-0.00001), -1)

    def test_modulo(self):
        """Test the modulo function."""
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
