#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Base class for tests."""

import unittest
import numpy as np
from numpy import angle
from tftb.utils import is_linear
from skimage.metrics import structural_similarity


class TestBase(unittest.TestCase):

    def assert_is_linear(self, signal, decimals=5):
        """Assert that the signal is linear."""
        self.assertTrue(is_linear(signal, decimals=decimals))

    def assert_is_analytic(self, signal, amlaw=None):
        """Assert that signal is analytic."""
        omega = angle(signal)
        if amlaw is not None:
            recons = np.exp(1j * omega) * amlaw
        else:
            recons = np.exp(1j * omega)
        real_identical = np.allclose(np.real(recons), np.real(signal))
        imag_identical = np.allclose(np.imag(recons), np.imag(signal))
        if not (imag_identical and real_identical):
            raise AssertionError("Signal is not analytic.")

    def assert_is_concave(self, signal):
        second_derivative = np.diff(np.diff(signal))
        if not np.all(second_derivative < 0):
            raise AssertionError("Signal is not concave.")

    def assert_is_convex(self, signal):
        second_derivative = np.diff(np.diff(signal))
        if not np.all(second_derivative > 0):
            raise AssertionError("Signal is not convex.")

    def assert_is_monotonic_increasing(self, signal):
        derivative = np.diff(signal)
        if not np.all(derivative >= 0):
            raise AssertionError("Signal is not monotonically increasing.")

    def assert_is_monotonic_decreasing(self, signal):
        derivative = np.diff(signal)
        if not np.all(derivative <= 0):
            raise AssertionError("Signal is not monotonically decreasing.")

    def assert_is_hermitian(self, x):
        """Assert that the input is a Hermitian matrix."""
        conj_trans = np.conj(x).T
        np.testing.assert_allclose(x, conj_trans)

    def assert_tfr_equal(self, x, y, sqmod=True, threshold=0.05, tol=0.99):
        """Assert that TFRs x and y are qualitatively equivalent."""
        if sqmod:
            x = np.abs(x) ** 2
            y = np.abs(y) ** 2
        x_thresh = np.amax(x) * threshold
        x[x <= x_thresh] = 0.0
        y_thresh = np.amax(y) * threshold
        y[y <= y_thresh] = 0.0
        x = np.ascontiguousarray(x)
        y = np.ascontiguousarray(y)
        similarity = structural_similarity(x, y, data_range=np.amax(x) - np.amin(x))
        self.assertTrue(similarity >= tol)
