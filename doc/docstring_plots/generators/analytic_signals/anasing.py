#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators.api import anasing
import matplotlib.pyplot as plt
import numpy as np

x = anasing(128)
plt.plot(np.real(x))
plt.grid()
plt.title('Lipschitz Singularity')
plt.show()
