#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.processing.freq_domain import inst_freq
from tftb.generators import fmlin


def plotifl(time_instants, iflaws, signal=None, **kwargs):
    """Plot normalized instantaneous frequency laws.

    :param time_instants: timestamps of the signal
    :param iflaws: instantaneous freqency law(s) of the signal.
    :param signal: if provided, display it.
    :type time_instants: array-like
    :type iflaws: array-like
    :type signal: array-like
    :return: None
    """
    indices = np.logical_not(np.isnan(iflaws))
    minif = np.amin(iflaws[indices])
    fig = plt.figure()
    if signal is not None:
        axsig = fig.add_subplot(211)
        axsig.set_position([0.10, 0.69, 0.80, 0.25])
        plt.sca(axsig)
        plt.plot(time_instants, np.real(signal))
        plt.title('Signal')
        plt.grid(True)
        plt.xlim(time_instants.min(), time_instants.max())
        axtfr = fig.add_subplot(212)
        axtfr.set_position([0.10, 0.21, 0.80, 0.45])
        plt.sca(axtfr)
    plt.plot(time_instants, iflaws)
    plt.xlim(time_instants.min(), time_instants.max())
    plt.grid(kwargs.get('grid', True))
    if minif >= 0:
        plt.ylim(0, 0.5)
    else:
        plt.ylim(-0.5, 0.5)
    plt.xlabel('Time')
    plt.ylabel('Normalized frequency')
    plt.title('Instantaneous frequency law(s)')
    plt.show()


if __name__ == '__main__':
    signal, _ = fmlin(256)
    time_samples = np.arange(3, 257)
    ifr = inst_freq(signal)
    plotifl(time_samples, ifr)
