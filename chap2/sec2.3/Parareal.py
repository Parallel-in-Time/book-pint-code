import numpy as np

def parareal(F, G, T, u0, N, K):
    """
    Implements the Parareal algorithm.

    Parameters
    ----------
    F : callable
        Fine solver, F(t0, t1, u0), advancing from t0 to t1 with initial condition u0.
    G : callable
        Coarse solver, G(t0, t1, u0), advancing from t0 to t1 with initial condition u0.
    T : float
        Final time.
    u0 : array-like
        Initial condition at t = 0.
    N : int
        Number of equidistant coarse time points.
    K : int
        Number of Parareal iterations.

    Returns
    -------
    U : list of numpy.ndarray
        List where U[k] contains the Parareal approximations at the coarse time points for each iteration k.
    """
    u0 = np.asarray(u0)
    TT = np.linspace(0, T, N+1)         # Coarse time mesh
    U = np.zeros((K+1, N+1, u0.size), dtype=u0.dtype)

    U[:, 0] = u0        # Initial condition (for all iterations)

    for n in range(N):  # Initial guess
        U[0, n+1] = G(TT[n], TT[n + 1], U[0, n])

    for k in range(K):  # Parareal iterations
        for n in range(N):
            U[k+1, n+1] = F(TT[n], TT[n+1], U[k, n]) \
                + G(TT[n], TT[n+1], U[k+1, n]) - G(TT[n], TT[n+1], U[k, n])

    return np.squeeze(U)
