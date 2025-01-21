import numpy as np

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
