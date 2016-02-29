#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

import unittest
from tftb.processing import time_domain as tmd
from tftb.generators import amplitude_modulated as am
from tftb.tests.base import TestBase


class TestTimeDomainProcessors(TestBase):

    def test_loctime(self):
        """Test computation of localized time characteristics."""
        signal = am.amgauss(160, 80, 50)
        tm, T = tmd.loctime(signal)
        self.assertAlmostEqual(tm, 79, places=6)
        self.assertAlmostEqual(T, 50, places=4)

if __name__ == '__main__':
    unittest.main()
