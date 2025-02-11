import numpy as np
from LagrangeWeights import lagrangeWeights

def PIDC(lam, tspan, u0, N, M, K):
    """
    Solves the Dahlquist equation with Parallel IDC.

    This function solves the Dahlquist equation with lambda on tspan=[tBeg, tEnd],
    starting from the initial value u0. It uses N parallel IDC time steps, with M points
    per time step (including the left and right boundary point). Then it computes the
    prediction (K=0) using Backward Euler, and corrects it using K Backward Euler sweeps.

    Parameters
    ----------
    lam : complex
        Parameter lambda in the Dahlquist equation.
    tspan : tuple
        Time interval (tBeg, tEnd).
    u0 : float
        Initial value.
    N : int
        Number of parallel IDC time steps.
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
    uStep = np.ones(K + 1, dtype=complex) * u0

    u = {}
    for n in range(N):
        u[n] = {}
        u[n][1] = np.zeros(M, dtype=complex)
        u[n][1][0] = uStep[0]

        # Initial guess using Backward Euler
        for m in range(M - 1):
            u[n][1][m + 1] = 1 / (1 - dt * lam) * u[n][1][m]

        # Sweep iterations
        for k in range(K):
            u[n][k + 2] = np.zeros(M, dtype=complex)
            u[n][k + 2][0] = uStep[k + 1]
            integrand = lam * u[n][k + 1]

            for m in range(M - 1):
                weights = lagrangeWeights(M)
                integral = (M - 1) * dt * sum(weights[m] * integrand)
                u[n][k + 2][m + 1] = 1 / (1 - dt * lam) * (
                    u[n][k + 2][m] - dt * lam * u[n][k + 1][m + 1] + integral
                )

        # Update initial values for the next step
        for k in range(K + 1):
            uStep[k] = u[n][k + 1][M - 1]

    # Combine each step
    uAll = [u0]
    for n in range(N):
        uAll.extend(u[n][K + 1][1:])

    return t, np.array(uAll)
