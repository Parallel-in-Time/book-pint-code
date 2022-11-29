# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 08:32:21 2018

@author: telu
"""
import numpy as np


def newton(f, fA, x0, nIter):
    """Find the zero of a function using the Newton method"""
    # Define variables
    eps = []
    xk = np.asarray(x0)

    # Compute first right hand side and error
    b = f(xk)
    eps.append(np.linalg.norm(b))

    # Newton loop
    for k in range(nIter):
        A = fA(xk)
        xkDiff = np.linalg.solve(A, b)
        xk -= xkDiff
        b = f(xk)
        eps.append(np.linalg.norm(b))

    return eps, xk
