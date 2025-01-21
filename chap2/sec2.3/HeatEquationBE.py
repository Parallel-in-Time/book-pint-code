import numpy as np
from scipy.sparse import spdiags, eye
from scipy.sparse.linalg import spsolve

def heatEquationBE(f, tspan, xspan, u0, gl, gr):
    """
    Solves the heat equation using the Backward Euler method.

    Parameters
    ----------
    f : callable
        Forcing function f(x, t).
    tspan : tuple
        Time interval (t0, T).
    xspan : tuple
        Spatial domain (x0, xL).
    u0 : ndarray
        Initial condition as a vector.
    gl : ndarray
        Dirichlet boundary condition on the left boundary.
    gr : ndarray
        Dirichlet boundary condition on the right boundary.

    Returns
    -------
    u : ndarray
        Solution matrix with dimensions (time steps, spatial points).
    """
    J = len(u0) - 1  # Number of spatial grid points minus 1
    N = len(gl) - 1  # Number of time steps minus 1
    dt = (tspan[1] - tspan[0]) / N
    dx = (xspan[1] - xspan[0]) / J
    x = np.linspace(xspan[0], xspan[1], J + 1)
    t = np.linspace(tspan[0], tspan[1], N + 1)
    
    u = np.zeros((N + 1, J + 1))  # Solution matrix
    u[0, :] = u0                 # Initial condition
    u[1:N + 1, 0] = gl[1:]       # Left boundary condition
    u[1:N + 1, J] = gr[1:]       # Right boundary condition

    # Construct the Laplacian matrix
    e = np.ones(J - 1)
    A = spdiags([e, -2 * e, e], [-1, 0, 1], J - 1, J - 1) / dx**2
    I = eye(J - 1)  # Identity matrix
    
    # Time-stepping loop
    for n in range(N):
        rhs = u[n, 1:-1] + dt * f(x[1:-1], t[n + 1])  # Source function
        rhs[0] += dt / dx**2 * gl[n + 1]              # Left boundary condition
        rhs[-1] += dt / dx**2 * gr[n + 1]             # Right boundary condition
        u[n + 1, 1:J] = spsolve(I - dt * A, rhs)      # Solve linear system
    
    return t, x, u
