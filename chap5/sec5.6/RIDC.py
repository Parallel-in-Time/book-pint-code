import numpy as np
from LagrangeWeights import lagrangeWeights
from IDC import IDC

def RIDC(lam, tspan, u0, N, M, K):
    """
    Solves the Dahlquist problem using RIDC.

    This function solves the Dahlquist equation with lambda on tspan=[tBeg, tEnd],
    starting from the initial value u0. It uses the equivalent of N IDC time steps,
    with M points per time step (including the left and right boundary point). Then it
    computes the prediction (K=0) using Backward Euler, and corrects using K Backward
    Euler sweeps.

    Parameters
    ----------
    lam : complex
        Parameter lambda in the Dahlquist equation.
    tspan : tuple
        Time interval (tBeg, tEnd).
    u0 : complex
        Initial value.
    N : int
        Number of IDC time steps.
    M : int
        Number of points per time step.
    K : int
        Number of correction sweeps.

    Returns
    -------
    t : array
        Time grid.
    u : array
        Solution values at time points.
    """
    nPoints = N * (M - 1) + 1
    t = np.linspace(tspan[0], tspan[1], nPoints)
    dt = t[1] - t[0]

    # Initialize solution array with complex dtype
    u = {1: np.zeros(nPoints, dtype=complex)}
    u[1][0] = u0

    # Initial guess using Backward Euler
    for n in range(nPoints - 1):
        u[1][n + 1] = 1 / (1 - dt * lam) * u[1][n]

    # IDC for the first M points
    for k in range(K):
        _, uI = IDC(lam, (t[0], t[M-1]), u0, 1, M, k)
        u[k + 2] = np.zeros(nPoints, dtype=complex)
        u[k + 2][:M] = uI

    # Loop on remaining points
    for n in range(M - 1, nPoints - 1):
        for k in range(K):
            integrand = lam * u[k + 1][n + 2 - M:n + 2]
            weights = lagrangeWeights(M)
            integral = (M - 1) * dt * np.sum(weights[-1] * integrand)
            u[k + 2][n + 1] = 1 / (1 - dt * lam) * (
                u[k + 2][n] - dt * lam * u[k + 1][n + 1] + integral
            )

    # Return only the last iterate
    return t, u[K + 1]
