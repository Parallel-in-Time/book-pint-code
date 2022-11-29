#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 15:16:38 2018

@author: lunet
"""
import numpy as np


def waveEquation(u0, ut0, c, L, tSpan, nStep):
    u0 = np.asarray(u0).ravel()
    nDOF = u0.size

    u = np.zeros((nDOF, nStep+1))
    t = np.linspace(0, tSpan, num=nStep+1)

    dx = L/(nDOF-1)
    dt = tSpan/nStep

    u[:, 0] = u0
    u[:, 1] = u0 + dt*ut0

    coeff = c*dt**2/dx**2
    for n in range(1, nStep):
        u[1:-1, n+1] = 2*u[1:-1, n] - u[1:-1, n-1] + coeff*(
            u[2:, n] - 2*u[1:-1, n] + u[0:-2, n])

    return t, u
