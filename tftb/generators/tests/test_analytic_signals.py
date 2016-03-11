#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for the analytic_signals module."""


import unittest
import numpy as np
from scipy import angle, unwrap
from scipy.stats import mode
from tftb.tests.base import TestBase
from tftb.generators import analytic_signals as ana


class TestAnalyticSignals(TestBase):

    def test_anaask(self):
        """Test analytic ASK signal."""
        signal, amlaw = ana.anaask(512, 64, 0.05)
        self.assert_is_analytic(signal, amlaw)
        np.testing.assert_allclose(np.abs(signal), amlaw)

    def test_anabpsk(self):
        """Test analytic BPSK signal."""
        signal, amlaw = ana.anabpsk(300, 30, 0.1)
        self.assertCountEqual(np.unique(amlaw), [-1, 1])
        self.assert_is_analytic(signal)

    def test_anafsk(self):
        """Test analytic phase shift keying."""
        signal, iflaw = ana.anafsk(512, 64, 3)
        self.assert_is_analytic(signal)

    def test_anapulse(self):
        """Test analytic unit impulse."""
        signal = ana.anapulse(512, 301)
        recons = np.zeros((512,))
        recons[301] = 1
        np.testing.assert_allclose(recons, np.real(signal), rtol=1e-5,
                                   atol=1e-5)

    def test_anaqpsk(self):
        """Test quaternary PSK signal."""
        signal, phases = ana.anaqpsk(512, 64, 0.25)
        self.assert_is_analytic(signal)
        # Count discontinuities in the signal and the phases and assert that
        # they appear in the same locations
        uphase = unwrap(angle(signal))
        dphase = np.diff(uphase)
        base_value = mode(dphase)[0][0]
        signal_phase_change = np.abs(dphase - base_value) > 0.0001
        ideal_phase_change = np.diff(phases) != 0
        np.testing.assert_allclose(signal_phase_change, ideal_phase_change)

    def test_anastep(self):
        """Test analytic unit step signal."""
        signal = ana.anastep(256, 128)
        recons = np.zeros((256,), dtype=float)
        recons[129:] = 1.0
        np.testing.assert_allclose(recons, np.real(signal), rtol=1e-5,
                                   atol=1e-5)

    def test_anasing(self):
        """Test the analytic singularity signal."""
        signal = ana.anasing(128)
        self.assert_is_analytic(signal, amlaw=np.abs(signal))


if __name__ == '__main__':
    unittest.main()
