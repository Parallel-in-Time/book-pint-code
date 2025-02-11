import numpy as np
import scipy.sparse as sp

def odeBEP(a, f, y0, t):
    """
    Solve ODE on a given time mesh with parallelized Backward Euler.

    This function solves the ordinary differential equation
    y' + a*y = f(t) with y(0) = y0 using Backward Euler on the time
    mesh t = [0, t1, ..., tN] by diagonalization, which could be performed in
    parallel. f here is a function.

    Parameters
    ----------
    a : float
        Coefficient of y in the ODE.
    f : callable
        Right-hand side function f(t).
    y0 : float
        Initial value y(0).
    t : array
        Time mesh.

    Returns
    -------
    y : array
        Solution values at time points t.
    """
    y = np.zeros(len(t))
    y[0] = y0
    N = len(t) - 1
    dt = np.diff(t)

    # Construct the matrix B
    B = sp.spdiags([-1 / np.append(dt[1:], 1), 1 / dt], [-1, 0], N, N).toarray()

    # Compute the eigen decomposition
    E, V = np.linalg.eig(B)

    # Evaluate f at the time points
    fv = f(t[1:])
    fv[0] += y0 / dt[0]

    # Solve the system using diagonalization
    E_diag = np.diag(E + a)
    y[1:] = (V @ np.linalg.solve(E_diag, np.linalg.solve(V, fv))).real

    return y
