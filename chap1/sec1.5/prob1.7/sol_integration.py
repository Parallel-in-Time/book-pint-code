#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 00:46:35 2018

@author: telu
"""
import numpy as np
import matplotlib.pyplot as plt
from solver import forwardEulerLin, backwardEulerLin

N = 64
T = 1/2
method = 'BE'

J = 99
h = 1/(J+1)
A = np.diag(np.ones(J-1), k=1) \
    + np.diag(np.ones(J-1), k=-1) \
    - 2*np.diag(np.ones(J), k=0)
A /= h**2

x = np.linspace(0, 1, num=J+2)[1:-1]

u0 = np.ones(J)
u0 *= 20

if method == 'FE':
    t, u = forwardEulerLin(A, u0, T, N)
elif method == 'BE':
    t, u = backwardEulerLin(A, u0, T, N)


fig = plt.figure('Solution')
plt.plot(x, u[:, N//4], '-', label='$T=1/8$', markevery=0.05)
plt.plot(x, u[:, N//2], '-', label='$T=1/4$', markevery=0.05)
plt.plot(x, u[:, N], '-', label='$T=1/2$', markevery=0.05)

nF = 5
pi = np.pi
n = np.arange(nF)
xTh = np.linspace(0, 1, num=500)
xMesh, nMesh = np.meshgrid(xTh, n)


def uTh(t):
    s = 80/pi/(2*nMesh+1) * np.exp(-(2*nMesh+1)**2 * pi**2 * t) \
        * np.sin((2*nMesh+1)*pi*xMesh)
    return np.sum(s, axis=0)


plt.plot(xTh, uTh(1/8), 'o', c='gray', markevery=0.05)
plt.plot(xTh, uTh(1/4), 'o', c='gray', markevery=0.05)
plt.plot(xTh, uTh(1/2), 'o', c='gray', markevery=0.05)

plt.xlabel('$x$')
plt.ylabel('$u(t)$')
plt.legend()
plt.grid(True)
fig.set_tight_layout(True)
plt.show()
