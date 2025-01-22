import numpy as np
from scipy.sparse import spdiags, eye, csr_matrix
from scipy.sparse.linalg import spsolve

def transportBE(f, a, tspan, xspan, u0, N):
    """
    Solves the periodic transport equation with Backward Euler.

    Parameters
    ----------
    f : callable
        Forcing function f(x, t).
    a : float
        Transport velocity.
    tspan : tuple
        Time interval (t0, T).
    xspan : tuple
        Spatial domain (x0, L).
    u0 : ndarray
        Initial condition as a vector.
    N : int
        Number of time steps.

    Returns
    -------
    u : ndarray
        Solution matrix with dimensions (time steps, spatial points).
    """
    J = len(u0) - 1  # Number of spatial grid points minus 1
    u = np.zeros((N + 1, J + 1))  # Solution matrix
    u[0, :] = u0  # Store initial condition

    dt = (tspan[1] - tspan[0]) / N  # Time step size
    dx = (xspan[1] - xspan[0]) / J  # Spatial step size
    x = np.linspace(xspan[0], xspan[1], J + 1)
    t = np.linspace(tspan[0], tspan[1], N + 1)

    # Construct the periodic centered first derivative matrix
    e = np.ones(J)
    A = spdiags([-e, e], [-1, 1], J, J) / (2 * dx)
    A = A.toarray()
    A[0, -1] = -1 / (2 * dx)  # Periodic boundary condition for the first row
    A[-1, 0] = 1 / (2 * dx)   # Periodic boundary condition for the last row

    I = eye(J)  # Identity matrix

    P = csr_matrix(I + a*dt*A)

    # Time-stepping loop
    for n in range(N):
        rhs = u[n, :-1] + dt*f(x[:-1], t[n + 1])  # Add source term
        u[n + 1, :-1] = spsolve(P, rhs)  # Solve linear system
        u[n + 1, -1] = u[n + 1, 0]  # Enforce periodic boundary condition

    return u
