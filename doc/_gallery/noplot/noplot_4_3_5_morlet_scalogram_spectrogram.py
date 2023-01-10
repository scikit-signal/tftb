#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Comparison of the Spectrogram and the morlet scalogram with their reassinged
counterparts.

Figure 4.36 from the tutorial
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import hamming

from tftb.generators import fmsin, fmlin, fmconst
from tftb.processing.cohen import Spectrogram
from tftb.processing.reassigned import spectrogram as re_spectrogram
from tftb.processing.reassigned import morlet_scalogram as re_morlet_scalogram
from tftb.processing import ideal_tfr


def generate_signal():
    sig1, if1 = fmsin(60, 0.16, 0.35, 50, 1, 0.35, 1)
    sig2, if2 = fmlin(60, 0.3, 0.1)
    sig3, if3 = fmconst(60, 0.4)

    sig = np.hstack((sig1, np.zeros((8,)), sig2 + sig3))
    iflaw = np.zeros((2, 128))
    iflaw[0, :] = np.hstack((if1, np.nan * np.ones((8,)), if2))
    iflaw[1, :] = np.hstack((np.nan * np.ones((68,)), if3))
    return sig, iflaw


def create_subplot(data, t=None, f=None, subplot=None, title=None,
                   xlabel=None, ylabel=None,
                   xtick_labels=None, ytick_labels=None,
                   kind="imshow"):

    plt.subplot(subplot)

    if kind == "contour":
        plt.contour(t, f, data, 1)
    elif kind == "imshow":
        plt.imshow(data, extent=[0, 128, 0, 0.5], aspect="auto", origin="lower")

    if xtick_labels is None:
        plt.gca().set_xticklabels([])
    if ytick_labels is None:
        plt.gca().set_yticklabels([])
    plt.grid(True)
    if title is not None:
        plt.title(title)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)


def example_4_3_5_morlet():
    sig, iflaw = generate_signal()
    tfr, t, f = ideal_tfr(iflaw)
    plt.figure(figsize=(10, 8))

    create_subplot(
        tfr,
        t=t,
        f=f,
        subplot=221,
        title="Ideal instantaneous frequencies",
        ylabel="Normalized Frequencies",
        ytick_labels=[],
        kind="contour")

    fs = 1.0
    nfft = 128
    nperseg = 33
    window = hamming(nperseg)
    noverlap = nperseg - 1
    return_onesided = False
    scaling = "density"
    mode = "psd"

    spec = Spectrogram(sig, n_fbins=nfft, fwindow=window)
    tfr, ts, freqs = spec.run(
        fs=fs,
        window=window,
        nperseg=nperseg,
        noverlap=noverlap,
        nfft=nfft,
        return_onesided=return_onesided,
        scaling=scaling,
        mode=mode
    )

    threshold = np.amax(np.abs(tfr)) * 0.05
    tfr[np.abs(tfr) <= threshold] = 0.0

    create_subplot(
        np.abs(tfr)[:64, :],
        subplot=222,
        title="Spectrogram",
    )

    _, tfr, _ = re_spectrogram(sig)
    tfr = tfr[:64, :]
    threshold = np.amax(np.abs(tfr) ** 2) * 0.05
    tfr[np.abs(tfr) ** 2 <= threshold] = 0.0

    create_subplot(
        np.abs(tfr) ** 2,
        subplot=223,
        title="Reassigned spectrogram",
        xlabel="Time",
        ylabel="Normalized Frequencies",
        xtick_labels=[],
        ytick_labels=[],
    )

    _, rtfr, _ = re_morlet_scalogram(sig)
    rtfr = rtfr[:64, :]
    threshold = np.amax(np.abs(rtfr) ** 2) * 0.05
    rtfr[np.abs(rtfr) ** 2 <= threshold] = 0.0

    create_subplot(
        np.abs(rtfr) ** 2,
        subplot=224,
        title="Reassigned Morlet Scalogram",
        xlabel="Time",
        xtick_labels=[],
    )

    plt.show()


if __name__ == "__main__":
    example_4_3_5_morlet()
