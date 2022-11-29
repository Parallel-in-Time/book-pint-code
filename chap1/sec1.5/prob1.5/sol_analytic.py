#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 01:10:08 2018

@author: telu
"""
import numpy as np
import matplotlib.pyplot as plt


nF = 5
pi = np.pi

n = np.arange(nF)
x = np.linspace(0, 1, num=10)
xMesh, nMesh = np.meshgrid(x, n)


# Function to compute the analytical solution at time t
def u(t):
    s = 80/pi/(2*nMesh+1) * np.exp(-(2*nMesh+1)**2 * pi**2 * t) \
        * np.sin((2*nMesh+1)*pi*xMesh)
    return np.sum(s, axis=0)


# Plot solution
plt.plot(x, u(1/8), label='$T=1/8$')
plt.plot(x, u(1/4), label='$T=1/8$')
plt.plot(x, u(1/2), label='$T=1/8$')
plt.legend()
plt.xlabel('$x$')
plt.ylabel('$u(x,T)$')
plt.grid(True)
plt.show()
