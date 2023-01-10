#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Tests for tftb.processing.cohen
"""

import unittest
import numpy as np
from scipy.signal.windows import kaiser
from tftb.processing import cohen
from tftb.generators import fmsin, fmlin
from tftb.tests.test_base import TestBase


class TestCohen(TestBase):

    def test_page_reality(self):
        """Test the reality property of the Page distribution."""
        signal, _ = fmsin(128)
        signal = signal / 128.0
        tfr, _, _ = cohen.PageRepresentation(signal).run()
        self.assertTrue(np.all(np.isreal(tfr)))

    def test_spectrogram_time_invariance(self):
        """Test the time invariance property of the spectrogram."""
        signal, _ = fmlin(128, 0.1, 0.4)
        window = kaiser(17, 3 * np.pi)
        tfr, ts, freqs = cohen.Spectrogram(signal, n_fbins=64, fwindow=window).run()
        shift = 64
        timeshifted_signal = np.roll(signal, shift)
        timeshifted_tfr, _, _ = cohen.Spectrogram(timeshifted_signal, n_fbins=64,
                                                  fwindow=window).run()
        rolled_tfr = np.roll(tfr, shift, axis=1)
        # the time invariance property holds mostly qualitatively. The shifted
        # TFR is not numerically indentical to the rolled TFR, having
        # differences at the edges; so clip with two TFRs where there are
        # discontinuities in the TFR.
        edge = 10
        xx = np.c_[timeshifted_tfr[:, edge:(shift - edge)],
                   timeshifted_tfr[:, (shift + edge):-edge]]
        yy = np.c_[rolled_tfr[:, edge:(shift - edge)],
                   rolled_tfr[:, (shift + edge):-edge]]
        np.testing.assert_allclose(xx, yy)

    def test_spectrogram_non_negativity(self):
        """Test that the spectrogram is non negative."""
        signal, _ = fmlin(128, 0.1, 0.4)
        window = kaiser(17, 3 * np.pi)
        tfr, _, _ = cohen.Spectrogram(signal, n_fbins=64, fwindow=window).run()
        self.assertTrue(np.all(tfr >= 0))

    def test_spectrogram_energy_conservation(self):
        """Test the energy conservation property of the spectrogram."""
        signal, _ = fmlin(128, 0.1, 0.4)
        window = kaiser(17, 3 * np.pi)
        tfr, ts, freqs = cohen.Spectrogram(signal, n_fbins=64, fwindow=window).run()
        e_sig = (np.abs(signal) ** 2).sum()
        self.assertAlmostEqual(tfr.sum().sum() / 64, e_sig)

    def test_spectrogram_reality(self):
        """Test the reality property of the spectrogram."""
        signal, _ = fmlin(128, 0.1, 0.4)
        window = kaiser(17, 3 * np.pi)
        tfr, _, _ = cohen.Spectrogram(signal, n_fbins=64, fwindow=window).run()
        self.assertTrue(np.all(np.isreal(tfr)))

    def test_spectrogram_linearity(self):
        """Test the linearity property of the spectrogram."""
        signal, _ = fmlin(128, 0.1, 0.4)
        window = kaiser(17, 3 * np.pi)
        tfr1, _, _ = cohen.Spectrogram(signal, n_fbins=64,
                                       fwindow=window).run()
        tfr2, _, _ = cohen.Spectrogram(signal * 2, n_fbins=64,
                                       fwindow=window).run()
        x = np.sum(np.sum(tfr2))
        y = np.sum(np.sum(tfr1))
        self.assertEqual(x / y, 4)

    def test_wigner_ville_energy(self):
        """Test the energy property of the Wigner Ville representation."""
        signal, _ = fmsin(128)
        signal = signal / 128.0
        tfr, _, _ = cohen.WignerVilleDistribution(signal).run()
        x = np.sum(np.sum(tfr))
        y = np.sum(np.abs(signal) ** 2) * 128
        self.assertEqual(x, y)

    def test_wigner_ville_projection(self):
        """Test the projection property of the Wigner Ville representation."""
        signal, _ = fmsin(128)
        tfr, _, _ = cohen.WignerVilleDistribution(signal).run()
        x = np.abs(signal) ** 2
        y = np.sum(tfr, axis=0) / 128
        np.testing.assert_allclose(x, y)

    def test_reality(self):
        """Test the reality property of the Wigner Ville representation."""
        signal, _ = fmsin(128)
        tfr, _, _ = cohen.WignerVilleDistribution(signal).run()
        self.assertTrue(np.all(np.isreal(tfr)))

    def test_wigner_ville_regionprops(self):
        """Test the regional property of the Wigner Ville representation."""
        signal, _ = fmsin(128)
        signal[64:] = 0
        tfr, _, _ = cohen.WignerVilleDistribution(signal).run()
        self.assertTrue(np.all(tfr[:, 64:] == 0))

        signal, _ = fmsin(128)
        signal[:64] = 0
        tfr, _, _ = cohen.WignerVilleDistribution(signal).run()
        self.assertTrue(np.all(tfr[:, :64] == 0))

    def test_pseudo_wv_energy(self):
        """Test the energy property of the pseudo WV representation."""
        signal, _ = fmsin(128)
        signal = signal / 128.0
        tfr, _, _ = cohen.PseudoWignerVilleDistribution(signal).run()
        x = np.sum(np.sum(tfr))
        y = np.sum(np.abs(signal) ** 2) * 128
        self.assertAlmostEqual(x, y, places=3)

    def test_pseudo_wv_reality(self):
        """Test the reality property of the pseudo WV representation."""
        signal, _ = fmsin(128)
        tfr, _, _ = cohen.PseudoWignerVilleDistribution(signal).run()
        self.assertTrue(np.all(np.isreal(tfr)))


if __name__ == '__main__':
    unittest.main()
