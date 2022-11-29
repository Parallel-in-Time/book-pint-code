#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 16:16:06 2018

@author: lunet
"""
import numpy as np
import matplotlib.pyplot as plt

# Main parameters
x0 = 2
nPas = 6
h = 1e-1
approx = True
approxType = 'FINDIFF'


# Function and derivative
def f(x):
    return np.sin(x)


def fPrim(x):
    if approx:
        if approxType == 'FINDIFF':
            return (f(x+h)-f(x-h))/(2*h)
        elif approxType == 'EVAL':
            return -0.9
    else:
        return np.cos(x)


# Newton with exact derivative
approx = False
xk = np.zeros(nPas+1)
xk[0] = x0
for k in range(nPas):
    xk[k+1] = xk[k] - f(xk[k])/fPrim(xk[k])
plt.semilogy(np.arange(nPas+1), abs(xk-np.pi), 'o-', label="f' exact")

# Newton with approximated constant derivative
approx = True
approxType = 'EVAL'
xk = np.zeros(nPas+1)
xk[0] = x0
for k in range(nPas):
    xk[k+1] = xk[k] - f(xk[k])/fPrim(xk[k])
plt.semilogy(np.arange(nPas+1), abs(xk-np.pi), 's-', label="f' const")

# Newton with finite difference derivative
approxType = 'FINDIFF'
xk = np.zeros(nPas+1)
xk[0] = x0
for k in range(nPas):
    xk[k+1] = xk[k] - f(xk[k])/fPrim(xk[k])
plt.semilogy(np.arange(nPas+1), abs(xk-np.pi), '^-', label="f' finDiff")

# Plotting parameters
plt.grid()
plt.xlabel('$N_{iter}$')
plt.ylabel('error')
plt.legend()
plt.savefig('newtonPi.pdf', bbox_inches='tight')
plt.tight_layout()
plt.show()
