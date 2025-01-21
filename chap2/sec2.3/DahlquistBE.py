import numpy as np

def dahlquistBE(la, tspan, u0, N):
    """
    Solves Dahlquist's test equation using the Backward Euler method.

    Parameters
    ----------
    la : float
        Parameter lambda in the test equation du/dt = la * u.
    tspan : tuple
        Time interval (start, end) for the solution.
    u0 : float
        Initial value of the solution.
    N : int
        Number of time steps.

    Returns
    -------
    t : numpy.ndarray
        Array of time points.
    u : numpy.ndarray
        Array of solution values at corresponding time points.
    """
    dt = (tspan[1] - tspan[0]) / N  # Time step size
    t = np.linspace(tspan[0], tspan[1], N + 1)  # Time points
    u = np.zeros(N + 1, dtype=complex)  # Solution array
    u[0] = u0  # Initial condition

    # Backward Euler method
    for n in range(N):
        u[n + 1] = u[n] / (1 - dt*la)

    return t, u
