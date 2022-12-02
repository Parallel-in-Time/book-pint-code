#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 15:31:38 2022

@author: telu
"""
import numpy as np
import matplotlib.pyplot as plt
from solver import forwardEuler, multipleShooting

# Lorentz parameters and functions
u0 = [20., 5., -5.]
sigma, rho, beta = 10., 28., 8/3

f = lambda t, u : [sigma*(u[1]-u[0]), u[0]*(rho-u[2])-u[1], u[0]*u[1]-beta*u[2]]

jac = lambda t, u : np.array([[-sigma, sigma, 0],
                              [rho-u[2], -1, -u[0]],
                              [u[1], u[0], -beta]])

# RHS of the coupled system (f and Jacobian)
def rhsFull(t, u):
    u, V = u[:3], u[3:].reshape((3, 3))
    uEval = f(t, u)
    jacEval = jac(t, u).dot(V)
    return np.ravel([uEval, *jacEval])

# Definition of propagator (solution and Jacobian) on subinterval
M = 10
def propagator(t0, t1, u0):
    u0 = np.ravel([u0, *np.eye(3)])  # Initial Jacobian is identity
    t, u = forwardEuler(rhsFull, t0, t1, u0, M)
    return u[-1, :3], u[-1, 3:].reshape((3, 3))  # rhs eval, Jacobian

# Multiple shooting parameters
N = 500
K = 9

# Utility plotting function
def setErrKPlot():
    plt.ylim(1e-13, 1e7)
    plt.ylabel(r'Error')
    plt.xlim(0, K)
    plt.xlabel(r'Iterations')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

# -----------------------------------------------------------------------------
# First experiment : Newton with different starting points
# -----------------------------------------------------------------------------

for tEnd in [0.5, 1, 1.5, 1.78, 2]:

    # Fine solution (for reference)
    tFine, uFine = forwardEuler(f, 0, tEnd, u0, N*M)
    uRef = uFine[::M]

    # Prediction for Multiple Shooting
    tPred, uPred = forwardEuler(f, 0, tEnd, u0, N)

    # Multiple Shooting (FE) solution
    times, uMS = multipleShooting(propagator, 0, tEnd, u0, N, K, uPred)

    # Plot error
    plt.figure('iterError')
    vecnorm = np.linalg.norm
    err = vecnorm(vecnorm(uMS-uRef, axis=-1), ord=np.inf, axis=-1)
    plt.semilogy(err, label=f'$T={tEnd}$')
    setErrKPlot()

    if tEnd in [2, 1]:

        plt.figure(f'timeError-{tEnd}')
        kMax = 6
        errNT = np.array(
            [np.linalg.norm(uk - uRef, axis=-1) for uk in uMS[:kMax]])
        for k, errk in enumerate(errNT):
            plt.semilogy(times, errk, label=f'K={k}')
        plt.legend()
        plt.ylim(1e-15, 1e4)
        plt.ylabel('Error')
        plt.xlabel('Time')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

plt.show()