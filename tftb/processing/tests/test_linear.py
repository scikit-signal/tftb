#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Tests for tftb.processing.linear
"""

from tftb.processing import linear
from tftb.generators import fmlin
from tftb.tests.base import TestBase
import numpy as np
import unittest


class TestLinear(TestBase):

    def test_stft_linearity(self):
        """Test the linearity property of the Fourier transform."""
        x = fmlin(128, 0.0, 0.2)[0]
        y = fmlin(128, 0.3, 0.5)[0]
        h = x + y
        tfr, _, _ = linear.ShortTimeFourierTransform(h).run()
        tfrx, _, _ = linear.ShortTimeFourierTransform(x).run()
        tfry, _, _ = linear.ShortTimeFourierTransform(y).run()
        np.testing.assert_allclose(tfr, tfrx + tfry)

    @unittest.skip("Known failure.")
    def test_stft_translation(self):
        """Test the time-shift property of the Fourier transform."""
        x = fmlin(128, 0.0, 0.2)[0]
        tfrx, _, freqs = linear.ShortTimeFourierTransform(x).run()
        x_shifted = np.roll(x, 64)
        tfrxs, _, _ = linear.ShortTimeFourierTransform(x_shifted).run()
        f_mul = np.exp(-1j * 2 * np.pi * 64 * freqs)
        np.testing.assert_allclose(tfrxs, f_mul * tfrx)

    def test_stft_modulation(self):
        """Test the modulation / frequency shifting property of STFT."""
        x = fmlin(128, 0.0, 0.2)[0]
        tfrx, _, _ = linear.ShortTimeFourierTransform(x).run()
        f_0 = 0.3
        h = np.exp(1j * 2 * np.pi * f_0 * np.arange(128)) * x
        tfrh, _, _ = linear.ShortTimeFourierTransform(h).run()
        tfrx_shifted = np.roll(tfrx, int(np.ceil(128 * 0.3)), axis=0)
        self.assert_tfr_equal(tfrx_shifted, tfrh, tol=0.95)

    @unittest.skip("Known failure")
    def test_stft_conjugation(self):
        x = fmlin(128, 0, 0.2)[0]
        h = np.conjugate(x)
        lhs, _, _ = linear.ShortTimeFourierTransform(h).run()
        rhs, _, _ = linear.ShortTimeFourierTransform(x[::-1]).run()
        rhs = np.conjugate(rhs)
        np.testing.assert_allclose(lhs, rhs)

if __name__ == '__main__':
    unittest.main()
