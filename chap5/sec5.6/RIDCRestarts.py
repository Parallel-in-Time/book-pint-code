import numpy as np
from RIDC import RIDC

def RIDCRestarts(lam, tspan, u0, N, M, K, R):
    """
    Solves the Dahlquist problem using RIDC with restarts.

    This function solves the Dahlquist test equation using RIDC with R restarts.

    Parameters
    ----------
    lam : complex
        Parameter lambda in the Dahlquist equation.
    tspan : tuple
        Time interval (tBeg, tEnd).
    u0 : complex
        Initial value.
    N : int
        Total number of IDC time steps.
    M : int
        Number of points per time step.
    K : int
        Number of correction sweeps.
    R : int
        Number of restarts.

    Returns
    -------
    t : array
        Time grid.
    u : array
        Solution values at time points.
    """
    NR = N // R
    t0 = tspan[0]
    dtR = (tspan[1] - tspan[0]) / R

    uR0 = u0
    t = np.array([t0])
    u = np.array([u0], dtype=complex)

    for r in range(R):
        tspanR = (t0 + r * dtR, t0 + (r + 1) * dtR)
        tR, uR = RIDC(lam, tspanR, uR0, NR, M, K)

        t = np.concatenate((t, tR[1:]))
        u = np.concatenate((u, uR[1:]))

        uR0 = u[-1]

    return t, u
