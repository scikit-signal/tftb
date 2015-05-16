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
from tftb.processing.utils import derive_window
import numpy as np
from scipy.signal import argrelmax, argrelmin, hanning


class TestMisc(TestBase):

    def test_window_derivative(self):
        window = hanning(210)
        derivative = derive_window(window)
        ix_win_maxima = np.argmax(window)
        self.assertAlmostEqual(derivative[ix_win_maxima], 0.0, places=3)

    def test_altes(self):
        ideal = np.array([0.00200822, -0.21928398, 0.66719239, 0.66719239,
                          -0.17666382, -0.17009953, -0.038399, -0.00083597])
        actual = misc.altes(8, 0.1, 0.5)
        np.testing.assert_allclose(ideal, actual, atol=1e-8, rtol=1e-8)

    def test_doppler(self):
        fm, am, iflaw = misc.doppler(512, 200.0, 65, 10, 50)
        self.assert_is_monotonic_decreasing(iflaw)

    def test_klauder(self):
        ideal = np.array([0.14899879, -0.16633309, -0.42806931, 0.16605633,
                          0.70769336, 0.16605633, -0.42806931, -0.16633309])
        actual = misc.klauder(8)
        np.testing.assert_allclose(ideal, actual, atol=1e-8, rtol=1e-8)

    def test_mexhat(self):
        ideal = np.array([-4.36444274e-09, -4.29488427e-04, -1.47862882e-01,
                          4.43113463e-01, -1.47862882e-01, -4.29488427e-04,
                          -4.36444274e-09])
        actual = misc.mexhat(0.5)
        np.testing.assert_allclose(ideal, actual, atol=1e-9, rtol=1e-9)
        maxima = argrelmax(actual)
        self.assertEqual(maxima[0].shape[0], 1)
        self.assertEqual(maxima[0][0], 3)
        minima = argrelmin(actual)
        self.assertEqual(minima[0].shape[0], 2)
        self.assertItemsEqual(minima[0], (2, 4))

    def test_gdpower(self):
        ideal_sig = np.array([0.11315600 + 0.j, 0.34703303 + 0.08691891j,
                              -0.02357698 + 0.49077882j, -0.34703303 + 0.03976496j,
                              -0.06600205 + 0.j, -0.34703303 - 0.03976496j,
                              -0.02357698 - 0.49077882j, 0.34703303 - 0.08691891j])
        ideal_f = np.array([0.125, 0.25, 0.375, 0.5])
        ideal_gpd = np.array([8.8, 6.45685425, 5.41880215, 4.8])
        ideals = (ideal_sig, ideal_gpd, ideal_f)
        actuals = misc.gdpower(8, 0.5)
        for i, ideal in enumerate(ideals):
            actual = actuals[i]
            np.testing.assert_allclose(ideal, actual, atol=1e-7, rtol=1e-7)


if __name__ == '__main__':
    unittest.main()
