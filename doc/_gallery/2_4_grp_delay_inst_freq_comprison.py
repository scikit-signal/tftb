#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 jaidev <jaidev@newton>
#
# Distributed under terms of the MIT license.

"""
Example in section 2.4 of the tutorial.
"""

import numpy as np
import matplotlib.pyplot as plt
from tftb.generators.api import amgauss, fmlin
from tftb.processing.api import loctime, locfreq, inst_freq, group_delay

"""
t=2:255;

sig1=amgauss(256,128,90).*fmlin(256,0,0.5);
[tm,T1]=loctime(sig1); [fm,B1]=locfreq(sig1);
T1*B1
---> T1*B1=15.9138
ifr1=instfreq(sig1,t); f1=linspace(0,0.5-1/256,256);
gd1=sgrpdlay(sig1,f1); plot(t,ifr1,’*’,gd1,f1,’-’)
sig2=amgauss(256,128,30).*fmlin(256,0.2,0.4);
[tm,T2]=loctime(sig2); [fm,B2]=locfreq(sig2);
T2*B2
---> T2*B2=1.224
ifr2=instfreq(sig2,t); f2=linspace(0.2,0.4,256);
gd2=sgrpdlay(sig2,f2); plot(t,ifr2,’*’,gd2,f2,’-’)
"""
time_instants = np.arange(2, 256)
sig1 = amgauss(256, 128, 90) * fmlin(256)[0]
tm, T1 = loctime(sig1)
fm, B1 = locfreq(sig1)
ifr1 = inst_freq(sig1, time_instants)
f1 = np.linspace(0, 0.5 - 1.0 / 256, 256)
gd1 = group_delay(sig1, f1)

plt.subplot(211)
plt.plot(time_instants, ifr1, '*', label='inst_freq')
plt.plot(gd1, f1, '-', label='group delay')
plt.xlim(0, 256)
plt.grid(True)
plt.legend()
plt.title("Time-Bandwidth product: {0}".format(T1 * B1))
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')


sig2 = amgauss(256, 128, 30) * fmlin(256, 0.2, 0.4)[0]
tm, T2 = loctime(sig2)
fm, B2 = locfreq(sig2)
ifr2 = inst_freq(sig2, time_instants)
f2 = np.linspace(0.02, 0.4, 256)
gd2 = group_delay(sig2, f2)


plt.subplot(212)
plt.plot(time_instants, ifr2, '*', label='inst_freq')
plt.plot(gd2, f2, '-', label='group delay')
plt.ylim(0.2, 0.4)
plt.xlim(0, 256)
plt.grid(True)
plt.legend()
plt.title("Time-Bandwidth product: {0}".format(T2 * B2))
plt.xlabel('Time')
plt.ylabel('Normalized Frequency')

plt.show()
