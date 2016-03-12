#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Tests for tftb.processing.postprocessing
"""

from tftb.tests.base import TestBase
from tftb.generators import atoms, fmlin
from tftb.processing import WignerVilleDistribution
from tftb.processing import postprocessing as pproc
from skimage.transform import hough_line_peaks, hough_line
import numpy as np
import unittest


class TestPostprocessing(TestBase):

    def test_renyi_information(self):
        """Check if Renyi entropy computation is correct."""
        sig = atoms(128, np.array([[64., 0.25, 20., 1.]]))
        tfr, _, _ = WignerVilleDistribution(sig).run()
        R1 = pproc.renyi_information(tfr)
        sig = atoms(128, np.array([[32., 0.25, 20., 1.],
                                   [96., 0.25, 20., 1.]]))
        tfr, _, _ = WignerVilleDistribution(sig).run()
        R2 = pproc.renyi_information(tfr)
        self.assertAlmostEqual(R2 - R1, 0.98, places=1)

    def test_ideal_tfr(self):
        """Test if the ideal TFR can be found using the instantaneous frequency
        laws."""
        _, iflaw1 = fmlin(128, 0.0, 0.2)
        _, iflaw2 = fmlin(128, 0.3, 0.5)
        iflaws = np.c_[iflaw1, iflaw2].T
        tfr, _, _ = pproc.ideal_tfr(iflaws)
        tfr[tfr == 1] = 255
        tfr = tfr.astype(np.uint8)
        hspace, angles, dists = hough_line(tfr)
        for x in hough_line_peaks(hspace, angles, dists):
            self.assertEqual(len(x), 2)


if __name__ == '__main__':
    unittest.main()
