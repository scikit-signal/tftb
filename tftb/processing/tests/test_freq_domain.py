#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""Tests for frequency domain processing functions."""


import unittest
from tftb.processing import freq_domain as fproc
from tftb.generators import frequency_modulated as fm
from tftb.tests.base import TestBase


class TestFrequencyDomainProcessing(TestBase):

    def test_instfreq(self):
        signal, _ = fm.fmlin(128, 0.05, 0.3, 50)
        ifreq = fproc.inst_freq(signal)
        self.assertAlmostEqual(ifreq.min(), 0.05, places=2)
        self.assertAlmostEqual(ifreq.max(), 0.3, places=2)
        self.assert_is_linear(ifreq)

    def test_locfreq(self):
        signal, _ = fm.fmlin(128, 0.05, 0.3, 50)
        input_avg_freq = (0.05 + 0.3) / 2.0
        avg_norm_freq, bandwidth = fproc.locfreq(signal)
        self.assertAlmostEqual(avg_norm_freq, input_avg_freq, places=2)


if __name__ == '__main__':
    unittest.main()
