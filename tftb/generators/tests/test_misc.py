#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for the generators.misc module."""

import unittest
from tftb.tests.test_base import TestBase
from tftb.generators import misc
from tftb.processing.utils import derive_window
import numpy as np
from scipy.signal import argrelmax, argrelmin
from scipy.signal.windows import hann


class TestMisc(TestBase):

    def test_window_derivative(self):
        """Test if the derivative of a window function is calculated
        properly."""
        window = hann(210)
        derivative = derive_window(window)
        ix_win_maxima = np.argmax(window)
        self.assertAlmostEqual(derivative[ix_win_maxima], 0.0, places=3)

    def test_altes(self):
        """Test the altes signal generation."""
        ideal = np.array([0.00200822, -0.21928398, 0.66719239, 0.66719239,
                          -0.17666382, -0.17009953, -0.038399, -0.00083597])
        actual = misc.altes(8, 0.1, 0.5)
        np.testing.assert_allclose(ideal, actual, atol=1e-8, rtol=1e-8)

    def test_doppler(self):
        """Test the doppler signal generation."""
        fm, am, iflaw = misc.doppler(512, 200.0, 65, 10, 50)
        self.assert_is_monotonic_decreasing(iflaw)

    def test_klauder(self):
        """Test the klauder wavelet generation."""
        ideal = np.array([0.14899879, -0.16633309, -0.42806931, 0.16605633,
                          0.70769336, 0.16605633, -0.42806931, -0.16633309])
        actual = misc.klauder(8)
        np.testing.assert_allclose(ideal, actual, atol=1e-8, rtol=1e-8)

    def test_mexhat(self):
        """Test the mexhat wavelet generation."""
        ideal = np.array([-4.36444274e-09, -4.29488427e-04, -1.47862882e-01,
                          4.43113463e-01, -1.47862882e-01, -4.29488427e-04,
                          -4.36444274e-09])
        actual = misc.mexhat(0.5)
        np.testing.assert_allclose(ideal, actual, atol=1e-9, rtol=1e-9)
        maxima = argrelmax(actual)
        self.assertEqual(maxima[0].shape[0], 1)
        self.assertEqual(maxima[0][0], 3)
        minima = argrelmin(actual)
        self.assertCountEqual(minima[0], (2, 4))

    def test_gdpower(self):
        """Test the gdpower generation."""
        ideal_sig = np.array([0.08540661 + 0.05077147j, 0.16735776 + 0.11542816j,
                              -0.08825763 + 0.17010894j, 0.04412953 - 0.01981114j,
                              -0.04981628 + 0.34985966j, -0.56798889 - 0.07983783j,
                              0.05266730 - 0.57074006j, 0.35650159 - 0.01577918j])
        ideal_f = np.array([0.125, 0.25, 0.375, 0.5])
        ideal_gpd = np.array([8.8, 6.45685425, 5.41880215, 4.8])
        ideals = (ideal_sig, ideal_gpd, ideal_f)
        actuals = misc.gdpower(len(ideal_sig), 0.5)
        for i, ideal in enumerate(ideals):
            actual = actuals[i]
            np.testing.assert_allclose(ideal, actual, atol=1e-7, rtol=1e-7)


if __name__ == '__main__':
    unittest.main()
