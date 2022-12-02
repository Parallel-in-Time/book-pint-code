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
    t0 : float
        The time of the initial solution
    tEnd : float
        The time of the final solution
    u0 : numpy vector or list,
        The initial solution
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


class LorenzSystem(object):
    r"""
    Instantiate the Lorenz system :

    .. math::
        \begin{align}
            \frac{dx}{dt} &= \sigma(y-x), \\
            \frac{dy}{dt} &= x(\rho-z)-y, \\
            \frac{dz}{dt} &=  xy-\beta z,
        \end{align}

    with :math:`\sigma`, :math:`\rho` and :math:`\beta` some given
    parameters.

    Parameters
    ----------
    sigma : float
        The :math:`\sigma` parameter (default=10)
    rho : float
        The :math:`\rho` parameter (default=28)
    beta : float
        The :math:`\beta` parameter (default=8/3)
    """
    def __init__(self, sigma=10, rho=28, beta=8/3):
        self.sigma = sigma
        self.rho = rho
        self.beta = beta

    def evalRHS(self, t, u):
        """Evaluate the right hand side of the Lorentz system"""
        x, y, z = u
        return self.sigma*(y-x), x*(self.rho-z)-y, x*y-self.beta*z

    def evalJacobian(self, t, u):
        """Evaluate the Jacobian of the Lorentz system"""
        x, y, z = u
        return np.array([[-self.sigma, self.sigma, 0],
                         [self.rho-z, -1, -x],
                         [y, x, -self.beta]
                        ])


def parareal(prop, rhs, t0, tEnd, u0, N, K, mF, mG):
    """
    Run the Parareal algorithm

    Parameters
    ----------
    prop : function(rhs, t0, tEnd, u0, nStep) -> t, u
        Time stepping function, that takes as argument :

        - rhs (function) : right hand side function (see below)
        - t0 and tEnd (float) : the integration interval bounds
        - u0 (vector or matrix) : the initial solution
        - nStep (int) : number of time integration steps

        It returns the time vector t of size (nStep+1)
        and the solution u as a (nStep+1, u0.shape) array.

    rhs : function(t, u) -> uEval
        The ODE right-hand-side function, that takes as argument :

        - t (float) : the time of evaluation
        - u (vector) : the solution for rhs evaluation

        It returns the solution uEval, same shape as u

    t0 : float
        Initial simulation time.
    tEnd : float
        End time of simulation.
    u0 : vector
        Initial solution of size nDOF.
    N : int
        Number of time sub-intervals.
    K : int
        Number of iterations.
    mF : int
        Number of time-steps for the fine solver per subintervals.
    mG : int
        Number of time-steps for the coarse solver per subintervals.

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
    u = np.zeros((K+1, N+1, nDOF), dtype=u0.dtype)

    # Time grid
    times = np.linspace(0, tEnd, N+1)

    # Set u^k_0 = u0
    u[:, 0] = u0

    # Prediction step using coarse solver
    tCoarse, uCoarse = prop(rhs, 0, tEnd, u0, N*mG)
    u[0, :] = uCoarse[::mG]

    # Iterations
    for k in range(K):
        for n in range(N):
            # Compute fine solution on subinterval
            t, uF = prop(rhs, times[n], times[n+1], u[k, n], mF)

            # Compute coarse solution from u_n^k
            t, uGk = prop(rhs, times[n], times[n+1], u[k, n], mG)

            # Compute coarse solution from u_n^{k+1}
            t, uGk1 = prop(rhs, times[n], times[n+1], u[k+1, n], mG)

            # Correction
            u[k+1, n+1] = uF[-1] + uGk1[-1] - uGk[-1]

    return times, u


def multipleShootingFE(rhs, jac, t0, tEnd, u0, N, K, M, uPred,
                       backTracking=False):
    """
    Run the MultipleShooting algorithm with "exact" Newton correction,
    using Forward-Euler time integration.

    Parameters
    ----------
    rhs : function(t, u) -> uEval
        The ODE right-hand-side function, that takes as argument :

        - t (float) : the time of evaluation
        - u (vector) : the solution for rhs evaluation

        It returns the solution uEval, same shape as u (nDOF)

    jac : function(t, u) -> Jf
        Function evaluating the Jacobian of the rhs, that takes as argument :

        - t (float) : the time of evaluation
        - u (vector) : the solution for Jacobian evaluation

        It returns the Jacobi matrix Jf, of shape (nDOF, nDOF)

    t0 : float
        Initial simulation time.
    tEnd : float
        End time of simulation.
    u0 : vector
        Initial solution of size nDOF.
    N : int
        Number of time sub-intervals.
    K : int
        Number of iterations.
    M : int
        Number of Forward-Euler time steps per subintervals.
    uPred : vector of size (nDOF) or matrix (N+1, nDOF)
        Prediction solution used for k=0.
    backTracking : bool, optional
        Wether or not use block adaptive back-tracking (experimental).
        The default is False.

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
    u = np.zeros((K+1, N+1, nDOF), dtype=u0.dtype)

    # Time grid
    times = np.linspace(0, tEnd, N+1)

    # Set u^k_0 = u0
    u[:, 0] = u0

    # Set u^0_n = uPred
    u[0, :] = uPred

    for k in range(K):
        for n in range(N):
            # Compute fine solution on subinterval
            t, uF = forwardEuler(rhs, times[n], times[n+1], u[k, n], M)

            # Store all inner fine solutions
            innerU = {tm: um for tm, um in zip(t, uF)}

            # Compute Jacobian
            def rhsJac(t, V):
                Jf = jac(t, innerU[t])
                return Jf.dot(V)

            t, V = forwardEuler(rhsJac, times[n], times[n+1], np.eye(nDOF), M)

            # Correction
            u[k+1, n+1] = uF[-1] + V[-1].dot(u[k+1, n] - u[k, n])


            if backTracking:
                Fk = u[k, n+1] - uF[-1]
                normK = np.linalg.norm(Fk)

                t, uFk1 = forwardEuler(
                    rhs, times[n], times[n+1], u[k+1, n], M)
                Fk1 = u[k+1, n+1] - uFk1[-1]

                a = 1.0
                beta = lambda a: a**(1/200)
                for i in range(20):
                    if np.linalg.norm(Fk1) < normK:
                        break
                    a /= 2
                    u[k+1, n+1] = a*uF[-1] \
                        + (1-a)*u[k, n+1] \
                        + beta(a)*V[-1].dot(u[k+1, n] - u[k, n])
                    Fk1 = u[k+1, n+1] - uFk1[-1]

    return times, u
