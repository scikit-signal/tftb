#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Comparison of the pseudo Page and the pseudo Margenau-Hill distributions with
their reassinged counterparts.

Figure 4.37 from the tutorial.
"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators import fmsin, fmlin, fmconst
from tftb.processing.cohen import PseudoPageRepresentation, PseudoMargenauHillDistribution
from tftb.processing.reassigned import pseudo_page as re_pseudo_page
from tftb.processing.reassigned import pseudo_margenau_hill as re_pseudo_margenau_hill

sig1, if1 = fmsin(60, 0.16, 0.35, 50, 1, 0.35, 1)
sig2, if2 = fmlin(60, 0.3, 0.1)
sig3, if3 = fmconst(60, 0.4)

sig = np.hstack((sig1, np.zeros((8,)), sig2 + sig3))
iflaw = np.zeros((2, 128))
iflaw[0, :] = np.hstack((if1, np.nan * np.ones((8,)), if2))
iflaw[1, :] = np.hstack((np.nan * np.ones((68,)), if3))

tfr, t, f = PseudoPageRepresentation(sig).run()
tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0

plt.figure(figsize=(10, 8))
plt.subplot(221)
plt.imshow(tfr[:64, :], extent=[0, 128, 0, 0.5], aspect='auto', origin='bottomleft')
plt.gca().set_xticklabels([])
plt.grid(True)
plt.title("Pseudo Page distribution")
plt.ylabel('Normalized Frequencies')

_, tfr, _ = re_pseudo_page(sig)
tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0
plt.subplot(222)
plt.imshow(tfr[:64, :], extent=[0, 128, 0, 0.5], aspect='auto', origin='bottomleft')
plt.grid(True)
plt.title("Reassigned Pseudo Page distribution")
plt.gca().set_xticklabels([])
plt.gca().set_yticklabels([])

tfr, _, _ = PseudoMargenauHillDistribution(sig).run()
tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0
plt.subplot(223)
plt.imshow(tfr[:64, :], extent=[0, 128, 0, 0.5], aspect='auto', origin='bottomleft')
plt.grid(True)
plt.title("Pseudo Margenau Hill distribution")
plt.xlabel('Time')
plt.ylabel('Normalized Frequencies')

_, rtfr, _ = re_pseudo_margenau_hill(sig)
tfr = np.abs(tfr) ** 2
threshold = np.amax(tfr) * 0.05
tfr[tfr <= threshold] = 0.0
plt.subplot(224)
plt.imshow(tfr[:64, :], extent=[0, 128, 0, 0.5], aspect='auto', origin='bottomleft')
plt.gca().set_yticklabels([])
plt.grid(True)
plt.title("Reassigned Pseudo Margenau Hill distribution")
plt.xlabel('Time')

plt.show()
