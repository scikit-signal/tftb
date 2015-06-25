#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Base time-frequency representation class.
"""

import numpy as np
import matplotlib.pyplot as plt


class BaseTFRepresentation(object):
    def __init__(self, signal, **kwargs):
        """Create a base time-frequency representation object.

        :param signal: Signal to be analyzed.
        :param **kwargs: Other arguments required for performing the analysis.
        :type signal: array-like
        :return: BaseTFRepresentation object
        :rtype:
        """
        if (signal.ndim == 2) and (1 in signal.shape):
            signal = signal.ravel()
        self.signal = signal
        timestamps = kwargs.get('timestamps')
        if timestamps is None:
            timestamps = np.arange(signal.shape[0])
        self.ts = self.timestamps = timestamps
        n_fbins = kwargs.get('n_fbins')
        if n_fbins is None:
            n_fbins = signal.shape[0]
        self.n_fbins = n_fbins
        fwindow = kwargs.get('fwindow')
        if fwindow is None:
            fwindow = self._make_window()
        self.fwindow = fwindow
        if self.n_fbins % 2 == 0:
            freqs = np.hstack((np.arange(self.n_fbins / 2),
                               np.arange(-self.n_fbins / 2, 0)))
        else:
            freqs = np.hstack((np.arange((self.n_fbins - 1) / 2),
                               np.arange(-(self.n_fbins - 1) / 2, 0)))
        self.freqs = freqs.astype(float) / self.n_fbins
        self.tfr = np.zeros((self.n_fbins, self.ts.shape[0]), dtype=complex)

    def _make_window(self):
        """Make a Hamming window function.

        The window function has a length equal to quarter of the length of the
        input signal.
        :return: Hamming window function.
        :rtype: array-like
        """

        h = np.floor(self.n_fbins / 4.0)
        h += 1 - np.remainder(h, 2)
        from scipy import hamming
        fwindow = hamming(int(h))
        fwindow = fwindow / np.linalg.norm(fwindow)
        return fwindow

    def plot(self, ax=None, kind='cmap', show=True, default_annotation=True,
             **kwargs):
        """Visualize the time frequency representation.

        :param ax: Axes object to draw the plot on.
        :param kind: One of "cmap" (default), "contour".
        :param show: Whether to call ``plt.show()``.
        :param default_annotation: Whether to make default annotations for the
            plot. Default annotations consist of setting the X and Y axis labels to
            "Time" and "Normalized Frequency" respectively, and setting the title
            to the name of the particular time-frequency distribution.
        :param **kwargs: Parameters to be passed to the plotting function.
        :type ax: matplotlib.axes.Axes object
        :type kind: str
        :type show: bool
        :type default_annotation: bool
        :return: None
        :rtype: None
        """
        extent = kwargs.pop('extent', None)
        if extent is None:
            extent = [self.ts.min(), self.ts.max(), self.freqs.min(),
                      self.freqs.max()]

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        if kind == "cmap":
            ax.imshow(self.tfr,
                      aspect='auto', origin='bottomleft', extent=extent,
                      **kwargs)
        else:
            t, f = np.meshgrid(self.ts, np.linspace(0, 0.5, self.signal.shape[0]))
            ax.contour(t, f, self.tfr, **kwargs)
        if default_annotation:
            ax.set_xlabel("Time")
            ax.set_ylabel("Normalized Frequency")
            ax.set_title(self.name.upper())
        if show:
            plt.show()
