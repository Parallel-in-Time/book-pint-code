import numpy as np

def nievergelt(solver, T, u0, N, uPred, Mn, width):
    """
    Implementation of Nievergelt's method for scalar ODEs.

    Parameters
    ----------
    solver : callable
        Function of the form solver(t0, t1, u0) that solves the ODE numerically 
        between t0 and t1 with u0 as the initial value, and returns the solution 
        over the interval [t0, t1].
    T : float
        The total time interval [0, T].
    u0 : float
        The initial value of the ODE.
    N : int
        The number of time sub-intervals dividing (0, T].
    uPred : ndarray
        Predicted solution, typically a cheap approximation.
    Mn : int
        The number of initial values for each sub-interval.
    width : float
        The width of the interval, centered around `U^0_n`, on which 
        the initial values are uniformly distributed.

    Returns
    -------
    U : ndarray
        The Nievergelt solution at the sub-interval endpoints.
    uTraj : ndarray
        The trajectories computed for each sub-interval.
    """
    dT = T / N
    TT = np.linspace(0, T, N + 1)
    uTraj = np.zeros((N, Mn, len(solver(TT[0], TT[1], u0))))

    # Compute first trajectory
    uTraj[0, 0, :] = solver(TT[0], TT[1], u0)

    # Compute remaining trajectories
    for n in range(1, N):
        uStart = uPred[n] + np.linspace(-width / 2, width / 2, Mn)
        for m in range(Mn):
            uTraj[n, m, :] = solver(TT[n], TT[n + 1], uStart[m])

    # Compute Nievergelt solution
    U = np.zeros(N + 1)
    U[0] = u0
    U[1] = uTraj[0, 0, -1]

    for n in range(1, N):
        m = findClosestInterval(U[n], uTraj[n, :, 0])
        uS1 = uTraj[n, m, 0]
        uS2 = uTraj[n, m + 1, 0]
        p = (U[n] - uS2) / (uS1 - uS2)
        U[n + 1] = p * uTraj[n, m, -1] + (1 - p) * uTraj[n, m + 1, -1]

    return U, uTraj


def findClosestInterval(x,y):
    """
    Finds the closest interval in a vector to a number.

    Parameters
    ----------
    x : float
        The value for which to find the closest interval.
    y : list or numpy.ndarray
        A sorted vector defining the intervals.

    Returns
    -------
    int
        The lower index `m` of the interval in `y` that contains `x`.
        If `x` is outside the range of `y`, returns the closest interval.
    """
    n = len(y)
    i = 0
    while i < n and x > y[i]:
        i += 1
    if i == 0:
        return 0
    elif i == n:
        return n-1
    else:
        return i-1