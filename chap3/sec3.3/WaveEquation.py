import numpy as np

def waveEquation(f, c, tspan, xspan, u0, u0t, gl, gr):
    """
    Solves the wave equation with centered finite differences.

    Parameters
    ----------
    f : callable
        Forcing function f(x, t).
    c : float
        Wave speed.
    tspan : list or tuple
        Time interval [t_start, t_end].
    xspan : list or tuple
        Spatial interval [x_start, x_end].
    u0 : ndarray
        Initial condition for u at t = t_start.
    u0t : ndarray
        Initial condition for du/dt at t = t_start.
    gl : ndarray
        Dirichlet boundary condition at the left boundary for all time steps.
    gr : ndarray
        Dirichlet boundary condition at the right boundary for all time steps.

    Returns
    -------
    u : ndarray
        Solution matrix with time along rows and space along columns.
    """
    J = len(u0) - 1  # Number of spatial grid points
    N = len(gl) - 1  # Number of time steps
    dt = (tspan[1] - tspan[0]) / N
    dx = (xspan[1] - xspan[0]) / J

    x = np.linspace(xspan[0], xspan[1], J + 1)
    t = np.linspace(tspan[0], tspan[1], N + 1)

    u = np.zeros((N + 1, J + 1))  # Solution matrix

    # Initial conditions
    u[0, :] = u0
    u[1, :] = u0 + dt * u0t  # Derivative by Taylor series expansion

    # Boundary conditions
    u[:, 0] = gl
    u[:, -1] = gr

    # Solve using centered finite differences
    for n in range(1, N):
        u[n + 1, 1:J] = (
            2 * u[n, 1:J] 
            - u[n - 1, 1:J]
            + (c * dt / dx) ** 2 * (u[n, 2:J + 1] - 2 * u[n, 1:J] + u[n, 0:J - 1])
            + dt ** 2 * f(x[1:-1], t[n])
        )

    return u
