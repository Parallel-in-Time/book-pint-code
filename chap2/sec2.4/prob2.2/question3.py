#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 15:31:38 2022

@author: telu
"""
import numpy as np
import matplotlib.pyplot as plt

from solver import LorenzSystem, forwardEuler, parareal, multipleShootingFE

# Initial parameters
u0 = [20., 5., -5.]
tEnd = 5
N = 500
M = 10

# Fine solution (for reference)
lorentz = LorenzSystem()
tFine, uFine = forwardEuler(lorentz.evalRHS, 0, tEnd, u0, N*M)
uRef = uFine[::M]

# Parareal solution
mG = 1
KP = 20
times, uP = parareal(forwardEuler, lorentz.evalRHS, 0, tEnd, u0, N, KP, M, mG)

# Utility plotting function
def setErrKPlot():
    plt.ylim(1e-13, 1e7)
    plt.ylabel(r'Error')
    plt.xlim(0, KP)
    plt.xlabel(r'Iterations')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

# -----------------------------------------------------------------------------
# First experiment : Parareal VS Newton
# -----------------------------------------------------------------------------

# Multiple Shooting (FE) solution
KN = 9
uPred = uP[0, :]
times, uN = multipleShootingFE(
    lorentz.evalRHS, lorentz.evalJacobian, 0, tEnd, u0, N, KN, M, uPred)

# Plot error
plt.figure()
errP = np.array(
    [np.linalg.norm(np.linalg.norm(uk - uRef, axis=-1),ord=np.inf, axis=-1)
     for uk in uP])
plt.semilogy(errP, label='Parareal error')
errN = np.array(
    [np.linalg.norm(np.linalg.norm(uk - uRef, axis=-1),ord=np.inf, axis=-1)
     for uk in uN])
plt.semilogy(errN, label='Multiple Shooting')
setErrKPlot()

plt.figure()
nK, kMax = 3, 10
errPT = np.array(
    [np.linalg.norm(uk - uRef, axis=-1) for uk in uP[:kMax:nK]])
for k, errk in enumerate(errPT):
    plt.semilogy(times, errk, '--', label=f'Parareal (K={k*nK})')
errNT = np.array(
    [np.linalg.norm(uk - uRef, axis=-1) for uk in uN[:kMax:nK]])
for k, errk in enumerate(errNT):
    plt.semilogy(times, errk, label=f'Multiple Shooting (K={k*nK})')
plt.legend()
plt.ylim(1e-13, 1e7)
plt.ylabel(r'Error')
plt.xlabel(r'Time')
plt.legend()
plt.grid(True)
plt.tight_layout()

# -----------------------------------------------------------------------------
# Second experiment : Multiple Shooting with different init. guesses
# -----------------------------------------------------------------------------

# Plot Parareal error
plt.figure()
errP = np.array(
    [np.linalg.norm(np.linalg.norm(uk - uRef, axis=-1), ord=np.inf, axis=-1)
     for uk in uP])
plt.semilogy(errP, label='Parareal')

# Multiple Shooting (FE) solution
for kS in [4, 6, 8]:
    KN = 8
    uPred = uP[kS, :]
    times, uN = multipleShootingFE(
        lorentz.evalRHS, lorentz.evalJacobian, 0, tEnd, u0, N, KN, M, uPred)

    # Plot Multiple Shooting error
    errN = np.array(
        [np.linalg.norm(np.linalg.norm(uk - uRef, axis=-1), ord=np.inf, axis=-1)
         for uk in uN])
    plt.semilogy(errN, label='Multiple Shooting, $K_{S}='f'{kS}$')
setErrKPlot()

# -----------------------------------------------------------------------------
# Third experiment : reducing simulation time (improve F and G accuracy)
# -----------------------------------------------------------------------------
tEnd = 1

# Fine solution (for reference)
lorentz = LorenzSystem()
tFine, uFine = forwardEuler(lorentz.evalRHS, 0, tEnd, u0, N*M)
uRef = uFine[::M]

# Parareal solution
times, uP = parareal(forwardEuler, lorentz.evalRHS, 0, tEnd, u0, N, KP, M, mG)

# Multiple Shooting (FE) solution
KN = 8
uPred = uP[0, :]
times, uN = multipleShootingFE(
    lorentz.evalRHS, lorentz.evalJacobian, 0, tEnd, u0, N, KN, M, uPred)

# Plot error
plt.figure()
errP = np.array(
    [np.linalg.norm(np.linalg.norm(uk - uRef, axis=-1),ord=np.inf, axis=-1)
     for uk in uP])
plt.semilogy(errP, label='Parareal error')
errN = np.array(
    [np.linalg.norm(np.linalg.norm(uk - uRef, axis=-1),ord=np.inf, axis=-1)
     for uk in uN])
plt.semilogy(errN, label='Multiple Shooting')
setErrKPlot()
