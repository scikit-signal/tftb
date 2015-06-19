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
            h = np.floor(n_fbins / 4.0)
            h += 1 - np.remainder(h, 2)
            from scipy import hamming
            fwindow = hamming(int(h))
        fwindow = fwindow / np.linalg.norm(fwindow)
        self.fwindow = fwindow
        self.tfr = np.zeros((self.n_fbins, self.ts.shape[0]), dtype=complex)

    def plot(self, ax=None, kind='cmap', **kwargs):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        if kind == "cmap":
            ax.imshow(self.tfr, extent=[self.ts.min(), self.ts.max(),
                                        self.freqs.min(), self.freqs.max()],
                      aspect='auto', origin='bottomleft', **kwargs)
        else:
            # FIXME: Implement contour plotting
            pass
