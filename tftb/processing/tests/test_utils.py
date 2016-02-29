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
Tests for tftb.processing.utils
"""

import unittest
import numpy as np
from tftb.processing import utils


class TestUtils(unittest.TestCase):

    def test_derive_window(self):
        from scipy.signal import gaussian
        g = gaussian(129, 10)
        dwindow = utils.derive_window(g)
        self.assertEqual(dwindow[64], 0)
        self.assertTrue(np.all(dwindow[:64] >= 0))
        self.assertTrue(np.all(dwindow[64:] <= 0))

if __name__ == '__main__':
    unittest.main()
