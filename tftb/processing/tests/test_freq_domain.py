#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for frequency domain processing functions."""


import unittest
import numpy as np
from tftb.processing import freq_domain as fproc
from tftb.generators import frequency_modulated as fm
from tftb.generators import amplitude_modulated as am
from tftb.tests.test_base import TestBase


class TestFrequencyDomainProcessing(TestBase):

    def test_instfreq(self):
        """Test instantaneous frequency calculation."""
        signal, _ = fm.fmlin(128, 0.05, 0.3, 50)
        ifreq = fproc.inst_freq(signal)[0]
        self.assertAlmostEqual(ifreq.min(), 0.05, places=2)
        self.assertAlmostEqual(ifreq.max(), 0.3, places=2)
        self.assert_is_linear(ifreq)

    def test_locfreq(self):
        """Test calculation of localized frequency characteristics."""
        signal, _ = fm.fmlin(128, 0.05, 0.3, 50)
        input_avg_freq = (0.05 + 0.3) / 2.0
        avg_norm_freq, bandwidth = fproc.locfreq(signal)
        self.assertAlmostEqual(avg_norm_freq, input_avg_freq, places=2)

    def test_group_delay(self):
        """Test Group delay calculation."""
        n_points = 128
        signal = am.amgauss(n_points, 64, 30) * fm.fmlin(n_points, 0.1, 0.4)[0]
        fnorm = np.arange(0.1, 0.38, step=0.04)
        grp_dlay = fproc.group_delay(signal, fnorm)
        self.assert_is_linear(grp_dlay, 0)


if __name__ == '__main__':
    unittest.main()
