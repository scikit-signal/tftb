#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Tests for tftb.processing.utils
"""

import unittest
import numpy as np
from tftb.processing import utils


class TestUtils(unittest.TestCase):

    def test_derive_window(self):
        """Test derivative of window function."""
        from scipy.signal import gaussian
        g = gaussian(129, 10)
        dwindow = utils.derive_window(g)
        self.assertEqual(dwindow[64], 0)
        self.assertTrue(np.all(dwindow[:64] >= 0))
        self.assertTrue(np.all(dwindow[64:] <= 0))

if __name__ == '__main__':
    unittest.main()
