import numpy as np

def forwardEuler(f, tspan, u0, N):
    """
    Solves a system of ODEs using the Forward Euler method.

    Parameters
    ----------
    f : callable
        Function defining the system of ODEs, f(t, u),
        has to return a np.ndarray
    tspan : tuple
        Time interval as (t0, t_end).
    u0 : array-like
        Initial condition.
    N : int
        Number of steps.

    Returns
    -------
    t : numpy.ndarray
        Array of time points.
    u : numpy.ndarray
        Solution array where each row corresponds to a time point.
    """
    dt = (tspan[1] - tspan[0]) / N  # Time step size
    t = np.linspace(tspan[0], tspan[1], N + 1)  # Time points
    u0 = np.asarray(u0)
    u = np.zeros((N + 1, len(u0)))  # Solution array
    u[0, :] = u0  # Set initial condition

    for n in range(N):  # Forward Euler steps
        u[n + 1, :] = u[n, :] + dt * f(t[n], u[n, :])

    return t, u
