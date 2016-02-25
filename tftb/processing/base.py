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
from scipy.interpolate import interp1d

TYPE1 = ["pseudo margenau-hill", "spectrogram", "pseudo page", "margenau-hill",
         "reassigned morlet scalogram", "gabor", "morlet scalogram",
         "reassigned pseudo margenau-hill", "reassigned spectrogram",
         "reassinged pseudo page", "reassigned gabor", "page", "rihaczek"]
TYPE2 = ["wigner-ville", "smoothed pseudo wigner-ville",
         "reassigned smoothed psedo wigner-ville", "d-flandrin", "bertrand",
         "unterberger", "pseudo wigner-ville", "reassigned psuedo wigner-ville",
         "scalogram"]
AFFINE = ["d-flandrin", "unterberger", "bertrand", "scalogram"]


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

    @property
    def has_negative_frequencies(self):
        return self.name.lower() in TYPE1

    @property
    def _isaffine(self):
        return self.name in AFFINE

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
        if extent is None:
            extent = [self.ts.min(), self.ts.max(), self.freqs.min(),
                      self.freqs.max()]
        if kind == "cmap":
            ax.imshow(self.tfr, cmap=cmap, origin="bottomleft", extent=extent,
                      aspect='auto')
        elif kind == "contour":
            if contour_x is None:
                contour_x = self.ts
            if contour_y is None:
                if show_tf:
                    if self.isaffine or self.name == "scalogram":
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

    def _plot_spectrum(self, ax, freqs, spec, scale, maxsp, freqr, nf2):
        # the axis here is the spectrum axis
        if scale == "linear":
            ax.plot(freqs, spec)
            ax.set_ylim(maxsp * self._viz_threshold, maxsp * 1.2)
        elif scale == "log":
            ax.plot(freqs, 10 * np.log10(spec / maxsp))
            ax.set_ylim(10 * np.log10(self._viz_threshold), 0)
        ax.set_xlim(self.fmin, self.fmax)
        if self._isaffine:
            freqs = np.linspace(freqr[0], freqr[nf2 - 1], nf2)
            spec = interp1d(0.5 * self.f_sample * np.arange(nf2) / nf2,
                            spec[:nf2])(freqs)
        else:
            freqs = freqr
            spec = spec[:nf2]
        maxsp = np.amax(spec)
        if scale == "linear":
            ax.set_ylim(maxsp * self._viz_threshold, maxsp * 1.2)
        elif scale == "log":
            ax.set_ylim(10 * np.log10(self._viz_threshold), 0)
        ax.set_xlim(self.fmin * self.f_sample, self.fmax * self.f_sample)

    def _annotate_spectrum(self, ax):
        ax.set_ylabel('Spectrum')
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.grid(True)
        if not self.isaffine:
            ax.invert_xaxis()
            ax.invert_yaxis()

    def plot(self, ax=None, kind='cmap', show=True, default_annotation=True,
             show_tf=False, scale="linear", threshold=0.05, f_sample=1.0,
             n_levels=64, fmin=0.0, fmax=None, **kwargs):
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
        self.f_sample = f_sample
        self.n_levels = n_levels
        self.fmin = fmin
        if fmax is None:
            fmax = 0.5 * self.f_sample
        self.fmax = fmax

        extent = kwargs.pop('extent', None)
        contour_x = kwargs.pop('contour_x', None)
        contour_y = kwargs.pop('contour_y', None)
        levels = kwargs.pop('levels', None)

        # data required for the plotting
        maxi = np.amax(self.tfr)
        tfrrow = self.tfr.shape[0]

        # Number of interesting frequency points to consider
        if self.has_negative_frequencies:
            nf2 = tfrrow / 2
        else:
            nf2 = tfrrow

        if self._isaffine:
            if "freq" not in kwargs:
                raise ValueError("Freq required for affine methods.")
            freq = kwargs.pop('freq')
        else:
            freq = 0.5 * np.arange(nf2) / nf2
        freqr = freq * f_sample
        self.ts = self.ts / f_sample

        if scale == "linear":
            if kind in ("surf", "mesh"):
                mini = np.amin(self.tfr)
            else:
                mini = np.max([np.amin(self.tfr), maxi * threshold])
            levels = np.linspace(mini, maxi, self.n_levels + 1)
        elif scale == "log":
            mini = np.max([np.amin(self.tfr), maxi * threshold])
            levels = np.logspace(np.log10(mini), np.log10(maxi), self.n_levels + 1)

        alpha = 2
        lt = self.ts.max() - self.ts.min() + 1
        if 2 * nf2 >= lt:
            spec = np.abs(np.fft.fft(self.signal[self.ts.min():self.ts.max()],
                                     2 * nf2)) ** 2
        else:
            nb_tranches_fog = np.floor(lt / (2 * nf2))
            spec = np.zeros((2 * nf2,))
            for i in range(int(nb_tranches_fog)):
                _spec = np.fft.fft(self.signal[self.ts.astype(int).min() + 2 * nf2 * i + np.arange(2 * nf2)])
                spec += np.abs(_spec) ** 2
            if lt > nb_tranches_fog * 2 * nf2:
                start = self.ts.min() + 2 * tfrrow * nb_tranches_fog
                spectre_fog = np.fft.fft(self.signal[start:self.ts.max()], 2 * nf2)
                spec += np.abs(spectre_fog) ** 2
        spec1 = np.sum(spec[:nf2])
        spec2 = np.sum(spec[nf2:])
        if spec2 > 0.1 * spec1:
            import warnings
            warnings.warn("The signal is not analytic.", UserWarning)
            if not np.all(np.isreal(self.signal)):
                alpha = 1

        if show_tf:
            if self._isaffine:
                f1 = freqr[0]
                f2 = freqr[nf2 - 1]
                d = f2 - f1
                nf4 = np.round((nf2 - 1) * self.f_sample / (2 * d)) + 1
                start = self.ts.astype(int).min()
                stop = self.ts.astype(int).max() + 1
                _spec_size = int(alpha * nf4)
                if _spec_size > spec.shape[0]:
                    spec = np.abs(np.fft.fft(self.signal[start:stop], int(alpha * nf4))) ** 2
                else:
                    spec[:_spec_size] = np.abs(np.fft.fft(self.signal[start:stop], int(alpha * nf4))) ** 2
                start = np.round(f1 * 2 * (nf4 - 1) / self.f_sample + 1)
                stop = np.round(f1 * 2 * (nf4 - 1) / self.f_sample + nf2) + 1
                spec = spec[start:stop]
                freqs = np.linspace(f1, f2, nf2)
            else:
                freqs = freqr
                spec = spec[:nf2]
            maxsp = np.amax(spec)

            fig, axTF = plt.subplots(figsize=(10, 8))
            self._plot_tfr(axTF, kind, extent, contour_x, contour_y, levels,
                        show_tf)
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(axTF)

            axTime = divider.append_axes("top", 1.2, pad=0.5)
            self._plot_signal(axTime)

            axSpec = divider.append_axes("left", 1.2, pad=0.5)
            self._plot_spectrum(axSpec, freqs, spec, scale, maxsp, freqr, nf2)

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
                from mpl_toolkits.mplot3d import Axes3D  # noqa
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

    def __tfrqview__(self, tfr=None, sig=None, t=None, method=None, **kwargs):
        if tfr is None:
            tfr = self.tfr
        tfrrow, tfrcol = tfr.shape
        if t is None:
            t = self.ts
        if method is None:
            method = "type1"
        if sig is None:
            sig = self.signal
        if self._isaffine:
            try:
                freq = kwargs['freq']
            except KeyError:
                raise KeyError("freq must be supplied to __tfrqview__ for affine reps.")
        else:
            freq = 0.5 * np.arange(nf2) / nf2  # noqa

        # Test of analyticity
        # FIXME: This test could be used for the unit testing too.
        lt_fog = t.max() - t.min() + 1
        nb_tranches_fog = np.floor(lt_fog / tfrrow)
        spec = np.zeros((tfrrow,))
        for i in range(nb_tranches_fog):
            _add = np.abs(np.fft.fft(sig[t.min() + tfrrow * i + np.arange(tfrrow)]))
            spec += _add ** 2
        if lt_fog > (nb_tranches_fog * tfrrow):
            sig_slice = sig[(t.min() + tfrrow * nb_tranches_fog):t.max()]
            spectre_fog = np.fft.fft(sig_slice, tfrrow)
            spec += np.abs(spectre_fog) ** 2
        spec1 = np.sum(spec[:(tfrrow / 2)])
        spec2 = np.sum(spec[(tfrrow / 2):])
        if spec2 > spec1 / 10.0:
            import warnings
            warnings.warn("The signal is not analytic", UserWarning)

        tfr = np.real(tfr)

    def __tfrview__(self, tfr=None, sig=None, t=None, method=None, kind="cmap",
            scale="linear", threshold=0.05, n_levels=64, nf2=None, fs=1.0,
            fmin=0.0, fmax=None, show_tf=True, **kwargs):
        if tfr is None:
            tfr = self.tfr
        tfrrow, tfrcol = tfr.shape
        if t is None:
            t = self.ts
        if method is None:
            method = "type1"
        if sig is None:
            sig = self.signal
        maxi = np.amax(tfr)

        # default params
        if nf2 is None:
            if self.has_negative_frequencies:
                nf2 = tfrrow / 2
            else:
                nf2 = tfrrow
        if fmax is None:
            fmax = fs * fmin

        # computation of isaffine and freq
        if self._isaffine:
            if "freq" not in kwargs:
                raise ValueError("Freq required for affine methods.")
            freq = kwargs.pop('freq')
        else:
            freq = 0.5 * np.arange(nf2) / nf2
        freqr = freq * fs
        ts = t / fs

        # update mini, levels, linlogstr, etc
        if scale == "linear":
            if kind in ("surf", "mesh"):
                mini = np.amin(tfr)
            else:
                mini = np.max([np.amin(tfr), maxi * threshold])
            levels = np.linspace(mini, maxi, n_levels + 1)
        elif scale == "log":
            mini = np.max([np.amin(tfr), maxi * threshold])
            levels = np.logspace(np.log10(mini), np.log10(maxi), n_levels + 1)

        # test of analyticity and computation of spec
        alpha = 2
        lt = t.max() - t.min() + 1
        if 2 * nf2 >= lt:
            spec = np.abs(np.fft.fft(sig[t.min():t.max()], 2 * nf2)) ** 2
        else:
            nb_tranches_fog = np.floot(lt / (2 * nf2))
            spec = np.zeros((2 * nf2,))
            for i in range(nb_tranches_fog):
                _spec = np.fft.fft(sig[t.min() + 2 * nf2 * i + np.arange(2 * nf2)])
                spec += np.abs(_spec) ** 2
            if lt > nb_tranches_fog * 2 * nf2:
                start = t.min() + 2 * tfrrow * nb_tranches_fog
                spectre_fog = np.fft.fft(sig[start:t.max()], 2 * nf2)
                spec += np.abs(spectre_fog) ** 2
        spec1 = np.sum(spec[:nf2])
        spec2 = np.sum(spec[nf2:])
        if spec2 > 0.1 * spec1:
            if not np.isreal(sig):
                alpha = 1

        if show_tf:
            if self._isaffine:
                f1 = freqr[0]
                f2 = freqr[nf2 - 1]
                d = f2 - f1
                nf4 = np.round((nf2 - 1) * fs / (2 * d)) + 1
                spec[:alpha * nf4] = np.abs(np.fft.fft(sig[t.min():t.max()],
                                                    alpha * nf4)) ** 2
                start = np.round(f1 * 2 * (nf4 - 1) / fs + 1)
                stop = np.round(f1 * 2 * (nf4 - 1) / fs + nf2)
                spec = spec[start:stop]
                freqs = np.linspace(f1, f2, nf2)
            else:
                freqs = freqr
                spec = spec[:nf2]
            maxsp = np.amax(spec)

            # the axis here is the spectrum axis
            if scale == "linear":
                plt.plot(freqs, spec)
                plt.ylim(maxsp * threshold * 0.01, maxsp * 1.2)
            elif scale == "log":
                plt.plot(freqs, 10 * np.log10(spec / maxsp))
                plt.ylim(10 * np.log10(threshold), 0)
            plt.xlim(fmin, fmax)
            if self._isaffine:
                freqs = np.linspace(freqr[0], freqr[nf2 - 1], nf2)
                spec = interp1d(0.5 * fs * np.arange(nf2) / nf2,
                                spec[:nf2])(freqs)
            else:
                freqs = freqr
                spec = spec[:nf2]
            maxsp = np.amax(spec)
            if scale == "linear":
                plt.ylim(maxsp * threshold, maxsp * 1.2)
            elif scale == "log":
                plt.ylim(10 * np.log10(threshold), 0)
            plt.xlim(fmin * fs, fmax * fs)

            # the axis here is the signal axis
            plt.plot(np.arange(t.min(), t.max()) / fs,
                     np.real(sig[t.min():t.max()]))
            plt.axis([t.min(), t.max(), np.amin(np.real(sig)),
                      np.amax(np.real(sig))])

            # The axis here is the TFR axis
            # for contour
            plt.contour(ts, freqr, tfr, levels)
            plt.ylim(fmin * fs, fmax * fs)

            # for images
            if scale == "linear":
                plt.imshow(ts, freqr, tfr)
            else:
                plt.imshow(ts, freqr, np.log10(tfr))
            plt.ylim(fmin * fs, fmax * fs)
