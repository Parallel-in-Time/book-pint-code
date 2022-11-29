#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 17:40:23 2018

@author: telu
"""
import numpy as np
import matplotlib.pyplot as plt

c = 1
x = np.linspace(-10, 10, num=500)


def u(x, t):
    return (np.exp(-(x+c*t)**2) + np.exp(-(x-c*t)**2))/2


plt.plot(x, u(x, 0), '--', label='$t=0$')
plt.plot(x, u(x, 1), 'o-', label='$t=1/c$', markevery=0.1)
plt.plot(x, u(x, 4), 's-', label='$t=4/c$', markevery=0.1)
plt.plot(x, u(x, 8), '^-', label='$t=8/c$', markevery=0.1)
plt.legend()
plt.xlabel('$x$')
plt.ylabel('$u(x,t)$')
plt.savefig('waveExp.pdf', bbox_inches='tight')
