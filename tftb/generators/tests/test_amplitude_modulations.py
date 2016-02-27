import unittest
import numpy as np
from numpy import pi
from scipy.signal import argrelmax
from tftb.tests.base import TestBase
import tftb.generators.amplitude_modulated as am


class TestAmplitudeModulated(TestBase):

    def test_amgauss(self):
        """Test if the gaussian amplitude modulator works correctly."""
        time_center = 63
        n_points = 128
        spread = 10
        signal = am.amgauss(n_points, time_center, spread)
        # parameters of the underlying gaussian function of the form
        # f(x) = a * exp( (-(x - b) **2) / (2 * (c ** 2)))
        a, b, c = 1, time_center, spread / np.sqrt(2 * pi)
        # Integral of a Gaussian is a * c * sqrt( 2 * pi)
        integral = a * c * np.sqrt(2 * pi)
        self.assertAlmostEqual(integral, signal.sum())

        # Other miscellaneous properties of a Gaussian
        maximum = argrelmax(signal)
        self.assertEqual(len(maximum), 1)
        self.assertEqual(maximum[0][0], time_center - 1)
        self.assertAlmostEqual(signal[time_center - 1], 1.0)

        self.assert_is_monotonic_increasing(signal[:(time_center - 1)])
        self.assert_is_monotonic_decreasing(signal[(time_center - 1):])

        infpl1 = np.floor(b - c).astype(int) - 1
        infpl2 = np.floor(b + c).astype(int)
        self.assert_is_convex(signal[:infpl1])
        self.assert_is_concave(signal[infpl1:infpl2])
        self.assert_is_convex(signal[infpl2:])

    def test_amexpos(self):
        """Test exponential amplitude modulation."""
        n_points, center, spread = 128, 63, 10.0
        one_sided = am.amexpos(n_points, center, spread, kind="unilateral")
        self.assertEqual(one_sided.max(), 1.0)
        self.assert_is_monotonic_decreasing(one_sided[center:])
        two_sided = am.amexpos(n_points, center, spread)
        self.assertEqual(two_sided.max(), 1.0)
        self.assert_is_monotonic_decreasing(two_sided[center:])
        self.assert_is_monotonic_increasing(two_sided[:center])

    def test_amrect(self):
        """Test rectangular amplitude modulation."""
        n_points, center, spread = 128, 63, 10.0
        signal = am.amrect(n_points, center, spread)
        self.assertEqual(signal.max(), 1.0)
        self.assertEqual(signal[center], 1.0)
        np.testing.assert_allclose(np.unique(signal), [0., 1.0])

    def test_amtriang(self):
        """Test triangular amplitude modulation."""
        n_points, center, spread = 128, 63, 10.0
        signal = am.amtriang(n_points, center, spread)
        self.assertEqual(signal.max(), 1.0)
        self.assertEqual(signal[center], 1.0)
        self.assert_is_monotonic_decreasing(signal[center:])
        self.assert_is_monotonic_increasing(signal[:center])


if __name__ == "__main__":
    unittest.main()
