#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 00:27:36 2018

@author: telu
"""
import numpy as np
import scipy.linalg as spl


def forwardEulerLin(A, u0, T, N):
    """
    Solve a linear system of ODE using Forward Euler.
    Considering the problem

    .. math::
        \\frac{dU}{dt} = AU,

    computes the numerical solution with Forward Euler, between :math:`t=0`
    and **T**, with **u0** the initial solution, using **N** steps.

    Parameters
    ----------
    A : matrix of size JxJ
        The matrix of the linear system
    u0 : vector of size J
        The initial vector
    T : float
        The final time of the solution
    N : int
        The number of numerical time steps

    Returns
    -------
    u : matrix of size JxN
        The solution at each time steps (including initial solution)
    t : vector of size N
        The times of the solutions
    """
    dt = T/N
    u0 = np.asarray(u0)
    J = u0.size
    u = np.zeros((J, N+1), dtype=u0.dtype, order='F')
    u[:, 0] = u0
    for j in range(N):
        np.dot(A, u[:, j], out=u[:, j+1])
        u[:, j+1] *= dt
        u[:, j+1] += u[:, j]
    t = np.linspace(0, T, N+1)
    return t, u


def backwardEulerLin(A, u0, T, N):
    """
    Solve a linear system of ODE using Backward Euler.
    Considering the problem

    .. math::
        \\frac{dU}{dt} = AU,

    computes the numerical solution with Backward Euler, between :math:`t=0`
    and **T**, with **u0** the initial solution, using **N** steps.

    Parameters
    ----------
    A : matrix of size JxJ
        The matrix of the linear system
    u0 : vector of size J
        The initial vector
    T : float
        The final time of the solution
    N : int
        The number of numerical time steps

    Returns
    -------
    u : matrix of size JxN
        The solution at each time steps (including initial solution)
    t : vector of size N
        The times of the solutions
    """
    dt = T/N
    u0 = np.asarray(u0)
    J = u0.size
    u = np.zeros((J, N+1), dtype=u0.dtype, order='F')
    u[:, 0] = u0
    R = np.copy(A)
    R *= -dt
    R += np.eye(J)
    for i in range(N):
        u[:, i+1] = np.linalg.solve(R, u[:, i])
    t = np.linspace(0, T, N+1)
    return t, u


def finDiffMatrixD2C2(J, L):
    """Compute the finite-difference matrix for the second derivative
    in space (D2) using centered finite differences of order 2.

    Parameters
    ----------
    J : int
        Number of mesh point in space (excluding boundary conditions)
    L : float
        Length of the domain (including boundary conditions)
    """
    h = L/(J+1)
    A = spl.toeplitz([-2., 1.]+(J-2)*[0.])
    A /= h**2
    return A
