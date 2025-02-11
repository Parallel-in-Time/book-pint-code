import numpy as np
from LagrangeWeights import lagrangeWeights

def IDC(lam, tspan, u0, N, M, K):
    """
    Solves the Dahlquist problem using Integral Deferred Correction.

    This function solves the Dahlquist equation with lambda on tspan=[tBeg, tEnd],
    starting from initial value u0. It uses N time steps for IDC, with M points per
    time step (including the left and right boundary point). Then it computes the
    prediction (K=0) using Backward Euler, and corrects using K Backward Euler sweeps.

    Parameters
    ----------
    lam : float
        Parameter lambda in the Dahlquist equation.
    tspan : tuple
        Time interval (tBeg, tEnd).
    u0 : float
        Initial value.
    N : int
        Number of time steps.
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
    uStep = u0

    u = {}
    for n in range(N):
        u[n] = {}
        u[n][1] = np.zeros(M, dtype=complex)
        u[n][1][0] = uStep

        # Initial guess using Backward Euler
        for m in range(M - 1):
            u[n][1][m + 1] = 1 / (1 - dt * lam) * u[n][1][m]

        # Iteration sweeps
        for k in range(K):
            u[n][k + 2] = np.zeros(M, dtype=complex)
            u[n][k + 2][0] = uStep
            integrand = lam * u[n][k + 1]

            for m in range(M - 1):
                weights = lagrangeWeights(M)
                integral = (M - 1) * dt * sum(weights[m] * integrand)
                u[n][k + 2][m + 1] = 1 / (1 - dt * lam) * (
                    u[n][k + 2][m] - dt * lam * u[n][k + 1][m + 1] + integral
                )

        uStep = u[n][K + 1][M - 1]

    # Combine each step solutions
    uAll = [u0]
    for n in range(N):
        uAll.extend(u[n][K + 1][1:])

    return t, np.array(uAll)
