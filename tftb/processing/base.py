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

    isaffine = False

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

    def _get_spectrum(self):
        if not self.isaffine:
            return np.fft.fftshift(np.abs(np.fft.fft(self.signal)) ** 2)
        nf2 = self.tfr.shape[0]
        spec = np.abs(np.fft.fft(self.signal[self.ts.min():(self.ts.max() + 1)],
                                 2 * nf2)) ** 2
        return spec[:nf2]

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
        # No need to normalize the window
        # fwindow = fwindow / np.linalg.norm(fwindow)
        return fwindow

    def _plot_tfr(self, ax, kind, extent, contour_x=None, contour_y=None,
                  levels=None, show_tf=True, cmap=plt.cm.gray):
        if kind == "cmap":
            ax.imshow(self.tfr, cmap=cmap, origin="bottomleft", extent=extent,
                      aspect='auto')
        elif kind == "contour":
            if contour_x is None:
                contour_x = self.ts
            if contour_y is None:
                if show_tf:
                    if self.isaffine:
                        contour_y = np.linspace(self.fmin, self.fmax, self.n_voices)
                    else:
                        contour_y = np.linspace(0, 0.5, self.signal.shape[0])
                else:
                    contour_y = np.linspace(0, 0.5, self.tfr.shape[0])
            contour_x, contour_y = np.meshgrid(contour_x, contour_y)
            if levels is not None:
                ax.contour(contour_x, contour_y, self.tfr, levels)
            else:
                if self.isaffine:
                    maxi = np.amax(self.tfr)
                    mini = max(np.amin(self.tfr), maxi * self._viz_threshold)
                    levels = np.linspace(mini, maxi, 65)
                    ax.contour(contour_x, contour_y, self.tfr, levels=levels)
                else:
                    ax.contour(contour_x, contour_y, self.tfr)

    def _annotate_tfr(self, ax):
        ax.grid(True)
        ax.set_xlabel("Time")
        ax.set_ylabel("Normalized Frequency")
        ax.set_title(self.name.upper())
        ax.yaxis.set_label_position("right")

    def _plot_signal(self, ax):
        ax.plot(np.real(self.signal))

    def _annotate_signal(self, ax):
        ax.set_xticklabels([])
        ax.set_xlim(0, self.signal.shape[0])
        ax.set_ylabel('Real part')
        ax.set_title('Signal in time')
        ax.grid(True)

    def _plot_spectrum(self, ax, freq_x, freq_y):
        k = int(np.floor(self.signal.shape[0] / 2.0))
        if freq_x is None:
            freq_x = self._get_spectrum()[::-1][:k]
        if freq_y is None:
            if self.isaffine:
                freq_y = self.freqs
            else:
                freq_y = np.arange(k)
        ax.plot(freq_x, freq_y)
        if not self.isaffine:
            ax.set_ylim(0, freq_y.shape[0] - 1)
        else:
            ax.set_ylim(freq_y[0], freq_y[-1])

    def _annotate_spectrum(self, ax):
        ax.set_ylabel('Spectrum')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.grid(True)
        if not self.isaffine:
            ax.invert_xaxis()
            ax.invert_yaxis()

    def plot(self, ax=None, kind='cmap', show=True, default_annotation=True,
             show_tf=False, scale="linear", threshold=0.05, **kwargs):
        """Visualize the time frequency representation.

        :param ax: Axes object to draw the plot on.
        :param kind: One of "cmap" (default), "contour".
        :param show: Whether to call ``plt.show()``.
        :param default_annotation: Whether to make default annotations for the
            plot. Default annotations consist of setting the X and Y axis labels to
            "Time" and "Normalized Frequency" respectively, and setting the title
            to the name of the particular time-frequency distribution.
        :param show_tf: Whether to show the signal and it's spectrum alongwith
            the plot. In this is True, the ``ax`` argument is ignored.
        :param **kwargs: Parameters to be passed to the plotting function.
        :type ax: matplotlib.axes.Axes object
        :type kind: str
        :type show: bool
        :type default_annotation: bool
        :return: None
        :rtype: None
        """
        self._viz_threshold = threshold

        extent = kwargs.pop('extent', None)
        if extent is None:
            extent = [self.ts.min(), self.ts.max(), self.freqs.min(),
                      self.freqs.max()]
        contour_x = kwargs.pop('contour_x', None)
        contour_y = kwargs.pop('contour_y', None)
        levels = kwargs.pop('levels', None)
        freq_x = kwargs.pop('freq_x', None)
        freq_y = kwargs.pop('freq_y', None)
        cmap = kwargs.pop("cmap", plt.cm.gray)

        if show_tf:
            fig, axTF = plt.subplots(figsize=(10, 8))
            self._plot_tfr(axTF, kind, extent, contour_x, contour_y, levels,
                           show_tf, cmap=cmap)
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(axTF)

            axTime = divider.append_axes("top", 1.2, pad=0.5)
            self._plot_signal(axTime)

            axSpec = divider.append_axes("left", 1.2, pad=0.5)
            self._plot_spectrum(axSpec, freq_x, freq_y)

            if default_annotation:
                self._annotate_tfr(axTF)
                self._annotate_signal(axTime)
                self._annotate_spectrum(axSpec)
        else:
            if (ax is None) and (kind != "surf"):
                fig = plt.figure()
                ax = fig.add_subplot(111)
            if kind == "cmap":
                ax.imshow(self.tfr,
                          aspect='auto', origin='bottomleft', extent=extent,
                          **kwargs)
            elif kind == "surf":
                from mpl_toolkits.mplot3d import Axes3D
                fig = plt.figure()
                ax = fig.gca(projection="3d")
                x = np.arange(self.signal.shape[0])
                y = np.linspace(0, 0.5, self.signal.shape[0])
                X, Y = np.meshgrid(x, y)
                ax.plot_surface(X, Y, np.abs(self.tfr), cmap=plt.cm.jet)
                if default_annotation:
                    ax.set_zlabel("Amplitude")
            elif kind == "wireframe":
                from mpl_toolkits.mplot3d import Axes3D  # NOQA
                ax = fig.gca(projection="3d")
                x = np.arange(self.signal.shape[0])
                y = np.linspace(0, 0.5, self.signal.shape[0])
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, np.abs(self.tfr), cmap=plt.cm.jet,
                                  rstride=3, cstride=3)
            else:
                t, f = np.meshgrid(self.ts, np.linspace(0, 0.5, self.tfr.shape[0]))
                ax.contour(t, f, self.tfr, **kwargs)
            if default_annotation:
                grid = kwargs.get('grid', True)
                ax.grid(grid)
                ax.set_xlabel("Time")
                ax.set_ylabel("Normalized Frequency")
                ax.set_title(self.name.upper())
        if show:
            plt.show()
