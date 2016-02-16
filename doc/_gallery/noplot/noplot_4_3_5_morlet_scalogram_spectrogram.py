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
from tftb.generators import fmsin, fmlin, fmconst
from tftb.processing.cohen import Spectrogram
from tftb.processing.reassigned import spectrogram as re_spectrogram
from tftb.processing.reassigned import morlet_scalogram as re_morlet_scalogram
from tftb.processing import ideal_tfr

sig1, if1 = fmsin(60, 0.16, 0.35, 50, 1, 0.35, 1)
sig2, if2 = fmlin(60, 0.3, 0.1)
sig3, if3 = fmconst(60, 0.4)

sig = np.hstack((sig1, np.zeros((8,)), sig2 + sig3))
iflaw = np.zeros((2, 128))
iflaw[0, :] = np.hstack((if1, np.nan * np.ones((8,)), if2))
iflaw[1, :] = np.hstack((np.nan * np.ones((68,)), if3))

tfr, t, f = ideal_tfr(iflaw)
plt.figure(figsize=(10, 8))
plt.subplot(221)
plt.contour(t, f, tfr, 1)
plt.grid(True)
plt.gca().set_xticklabels([])
plt.title("Ideal instantaneous frequencies")
plt.ylabel('Normalized Frequencies')

tfr, _, _ = Spectrogram(sig).run()
threshold = np.amax(np.abs(tfr)) * 0.05
tfr[np.abs(tfr) <= threshold] = 0.0
plt.subplot(222)
plt.imshow(np.abs(tfr)[:64, :], extent=[0, 128, 0, 0.5], aspect='auto', origin='bottomleft')
plt.grid(True)
plt.gca().set_xticklabels([])
plt.gca().set_yticklabels([])
plt.title("Spectrogram")

_, tfr, _ = re_spectrogram(sig)
tfr = tfr[:64, :]
threshold = np.amax(np.abs(tfr) ** 2) * 0.05
tfr[np.abs(tfr) ** 2 <= threshold] = 0.0
plt.subplot(223)
plt.imshow(np.abs(tfr) ** 2, extent=[0, 128, 0, 0.5], aspect='auto', origin='bottomleft')
plt.grid(True)
plt.title("Reassigned spectrogram")
plt.xlabel('Time')
plt.ylabel('Normalized Frequencies')

_, rtfr, _ = re_morlet_scalogram(sig)
rtfr = rtfr[:64, :]
threshold = np.amax(np.abs(rtfr) ** 2) * 0.05
rtfr[np.abs(rtfr) ** 2 <= threshold] = 0.0
plt.subplot(224)
plt.imshow(np.abs(rtfr) ** 2, extent=[0, 128, 0, 0.5], aspect='auto', origin='bottomleft')
plt.gca().set_yticklabels([])
plt.grid(True)
plt.title("Reassigned Morlet Scalogram")
plt.xlabel('Time')

plt.show()
