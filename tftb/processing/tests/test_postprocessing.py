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
Tests for tftb.processing.postprocessing
"""

from tftb.tests.base import TestBase
from tftb.generators import atoms
from tftb.processing import WignerVilleDistribution
from tftb.processing import postprocessing as pproc
import numpy as np
import unittest


class TestPostprocessing(TestBase):

    def test_renyi_information(self):
        sig = atoms(128, np.array([[64., 0.25, 20., 1.]]))
        tfr, _, _ = WignerVilleDistribution(sig).run()
        R1 = pproc.renyi_information(tfr)
        sig = atoms(128, np.array([[32., 0.25, 20., 1.],
                                   [96., 0.25, 20., 1.]]))
        tfr, _, _ = WignerVilleDistribution(sig).run()
        R2 = pproc.renyi_information(tfr)
        self.assertAlmostEqual(R2 - R1, 0.98, places=1)


if __name__ == '__main__':
    unittest.main()
