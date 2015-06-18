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
from tftb.generators.misc import atoms

coordinates = np.array([[32.0, 0.3, 32.0, 1.0],
                        [56.0, 0.15, 48.0, 1.22]])
sig = atoms(128, coordinates)
plt.plot(np.real(sig))
plt.grid()
plt.xlim(xmax=128)
plt.title('Gaussian Atoms')
plt.show()
