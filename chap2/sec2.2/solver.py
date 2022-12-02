#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 15:13:34 2022

@author: telu
"""
import numpy as np


def forwardEuler(f, t0, tEnd, u0, nStep):
    """
    Uses the Forward Euler method to solve the system of first order
    ordinary differential equations

    .. math::
        \\frac{du}{dt} = f(t,u)

    with :math:`u` a vector of size :math:`N_{dof}`.
    The final solution is computed using :math:`N_{step}` time steps.
    At each time steps, the updated solution is computed using the formula :

    .. math::
        u(t+dt) = u(t) + dt \\times f(t,u(t)),

    where :math:`dt=\\frac{t_{end}-t_0}{N_{step}}`.

    Parameters
    ----------
    f : function
        The :math:`f` operator, as a function that take a scalar **t** and
        a numpy vector **u** as argument,
        and returns a vector of the same size as **u**.
    t0, tEnd : float
        Initial and end simulation time.
    u0 : vector of size (nDOF)
        Initial solution of size nDOF.
    nStep : int
        The number of time step to be performed

    Returns
    -------
    t : numpy vector of size :math:`N_{step}+1`
        The discrete times from **t0** to **tEnd** when the solution was
        computed (including the initial time).
    u : numpy array of shape :math:`(N_{step}+1, u0.shape)`
        The solution of the ODE at each time steps, including the initial time.
    """
    # Create output variables
    u0 = np.asarray(u0)
    u = np.zeros((nStep+1,) + u0.shape)
    t = np.linspace(t0, tEnd, nStep+1)

    # Store initial solution and time
    u[0] = u0

    # Compute time-step
    dt = (tEnd-t0)/nStep

    # Loop on every time-step
    for i in range(nStep):
        # Evaluation of the operator
        u[i+1] = f(t[i], u[i])
        # Multiplication by time-step
        u[i+1] *= dt
        # Addition of previous step solution
        u[i+1] += u[i]

    return t, u


def multipleShooting(prop, t0, tEnd, u0, N, K, uPred):
    """
    Implementation of Multiple Shooting using "exact" Newton correction

    Parameters
    ----------
    prop : function(t0, t1, u0) -> u1, Jac1
        Time-propagator function for u0 between t0 and t1, that also
        propagates the Jacobian of the propagator with respect to u0
        evaluated on (t1, u1).
    t0, tEnd : float
        Initial and end simulation time.
    u0 : vector of size (nDOF)
        Initial solution of size nDOF.
    N : int
        Number of time sub-intervals.
    K : int
        Number of iterations.
    uPred : vector of size (nDOF) or matrix (N+1, nDOF)
        Prediction solution used for k=0.

    Returns
    -------
    times : vector of size (N+1)
        Sub-interval interface times.
    u : array of size (K+1, N+1, nDOF)
        Parareal solution for each iterations and time-subintervals.
    """
    # Initializing storage variables
    u0 = np.asarray(u0)
    nDOF = u0.size
    u = np.zeros((K+1, N+1, nDOF))

    # Time grid
    times = np.linspace(t0, tEnd, N+1)

    # Set
    u[:, 0] = u0

    # Set u^0_n = uPred
    u[0, :] = uPred

    for k in range(K):
        for n in range(N):
            # Compute fine solution and Jacobian
            uF, V = prop(times[n], times[n+1], u[k, n])
            # Correction update
            u[k+1, n+1] = uF + V.dot(u[k+1, n] - u[k, n])

    return times, u
