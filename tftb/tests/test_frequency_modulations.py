#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for the frequency_modulated module."""

import unittest
import numpy as np
from tftb.tests.base import TestBase
from tftb.generators import frequency_modulated as fm


class TestFrequencyModulated(TestBase):
    """Tests for the frequency_modulated module."""

    def test_fmconst(self):
        """Test constant frequency modulation."""
        n_points, fnorm, center = 129, 0.5, 64
        signal, iflaw = fm.fmconst(n_points, fnorm, center)
        real = np.real(signal)
        np.testing.assert_allclose(iflaw.ravel(), fnorm * np.ones((n_points,)))
        omega = np.arccos(real)[1:] / (2 * np.pi) / np.arange(n_points)[1:]
        self.assertEqual(omega.max(), fnorm)
        self.assert_is_analytic(signal)

    def test_fmhyp(self):
        """Test hyperbolic frequency modulation."""
        n_points, p1, p2 = 128, (1, 0.5), (32, 0.1)
        signal, iflaw = fm.fmhyp(n_points, p1, p2)
        self.assertEqual(iflaw[p2[0] - 1], p2[1])
        self.assert_is_analytic(signal)

    def test_fmlin(self):
        """Test linear frequency modulation."""
        n_points, init_freq, final_freq, center = 128, 0.05, 0.3, 50
        signal, iflaw = fm.fmlin(n_points, init_freq, final_freq, center)
        np.testing.assert_allclose(iflaw,
                                   np.linspace(init_freq, final_freq,
                                               n_points))
        self.assert_is_analytic(signal)

    def test_fmodany(self):
        """Test arbitrary modulation."""
        iflaw = np.pad(np.ones((32,)), pad_width=(48, 48), mode='constant',
                       constant_values=0)
        signal = fm.fmodany(iflaw / 2.0)
        self.assert_is_analytic(signal)
        np.testing.assert_allclose(np.real(signal)[:48],
                                   np.ones((48,)))
        np.testing.assert_allclose(np.real(signal)[-48:],
                                   np.ones((48,)))
        x = fm.fmconst(32, 0.5, 1)[0]
        np.testing.assert_allclose(np.real(x), np.real(signal)[48:48 + 32])

    def test_fmpar(self):
        """Test parabolic frequency modulation."""
        coeffs = (0.4, -0.0112, 8.6806e-05)
        n_points = 128
        signal, iflaw = fm.fmpar(n_points, coeffs)
        a, b, c = coeffs
        xx = np.arange(n_points)
        parabola = a + b * xx + c * (xx ** 2)
        np.testing.assert_allclose(parabola, iflaw, rtol=1e-3, atol=1e-3)
        self.assert_is_analytic(signal)

    def test_fmpower(self):
        """Test power law frequency modulation."""
        n_points = 128
        degree = 0.5
        coefficients = 1, 0.5, 100, 0.1
        signal, iflaw = fm.fmpower(n_points, degree, coefficients)
        self.assertAlmostEqual(iflaw[coefficients[2]], coefficients[3],
                               places=1)
        self.assert_is_analytic(signal)

    def test_fmsin(self):
        """Test sinusoidal frequency modulation."""
        n_points = 140
        min_freq, max_freq = 0.05, 0.45
        period = 100
        center = 20
        center_freq = 0.3
        freq_dir = -1
        signal, iflaw = fm.fmsin(n_points, min_freq, max_freq, period, center,
                                 center_freq, freq_dir)
        self.assert_is_analytic(signal)
        xx = np.arange(n_points) - center - 4
        f_sample = 1.0 / period
        yy = freq_dir * np.sin(2 * np.pi * f_sample * xx)
        yy = ((max_freq - min_freq) / 2) * yy + (max_freq + min_freq) / 2
        np.testing.assert_allclose(yy, iflaw, rtol=1e-3, atol=1e-3)

if __name__ == '__main__':
    unittest.main()
