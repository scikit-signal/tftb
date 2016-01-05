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
from scipy import angle
from tftb.utils import is_linear

# yoder:
# let's add at least some backwards python2.x compatibility for now.
import sys
#py_ver = sys.version_info.major
# and assume only backwards revisions (for now):
ispy2 = (sys.version_info.major<3)		# maybe we should work this into TestBase?

class TestBase(unittest.TestCase):

    # yoder: add __init__()
    def __init__(self, *args, **kwargs):
    	# handle various bits, including some python2-3 compatibility, then execute base __init__ as super()
    	if not ispy2:
    		# re-map some function calls:
    		self.assertItemsEqual = self.assertCountEqual		#(inherited from unittest.TestCase)
    		# ... and others...
    		#
    	# and let's go backwards as well, just to be sure (now, we can correct all the downstream code and remove this hack at a later time...).
    	if ispy2:
    		self.assertCountEqual = self.assertItemsEqual
    	#
    	super(TestBase,self).__init__(*args, **kwargs)
    
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
