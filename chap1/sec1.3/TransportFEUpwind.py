import numpy as np

def transportFEUpwind(f, u0, ua, a, b, J, T, N):
    """
    Solves the 1D transport equation u_t + u_x = f using upwind finite differences
    in space and forward Euler in time.

    Parameters
    ----------
    f : function
        Right-hand side function f(x, t).
    u0 : function
        Initial condition function u0(x).
    ua : function
        Boundary condition function ua(t) at x = a.
    a : float
        Left boundary of the spatial domain.
    b : float
        Right boundary of the spatial domain.
    J : int
        Number of subintervals in the spatial domain.
    T : float
        Final time of the simulation.
    N : int
        Number of time steps.

    Returns
    -------
    u : ndarray
        Solution array of shape (J+1, N+1).
    x : ndarray
        Spatial grid points.
    t : ndarray
        Time grid points.
    """
    dx = (b - a) / J
    dt = T / N
    x = np.linspace(a, b, J + 1)
    t = np.linspace(0, T, N + 1)
    u = np.zeros((J + 1, N + 1))
    
    # Apply boundary and initial conditions
    u[0, :] = ua(t)  # Boundary condition at x = a
    u[:, 0] = u0(x)  # Initial condition
    
    # Time-stepping loop
    for n in range(N):
        u[1:J+1, n+1] = (
            u[1:J+1, n]
            - dt / dx * (u[1:J+1, n] - u[0:J, n])
            + dt * f(x[1:J+1], t[n])
        )
    
    return u, x, t
