#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""

"""


from tftb.generators import anastep
import matplotlib.pyplot as plt
import numpy as np

x = anastep(256, 128)
plt.plot(np.real(x))
plt.xlim(0, 256)
plt.grid()
plt.title('Analytic Step Signal')
plt.show()
