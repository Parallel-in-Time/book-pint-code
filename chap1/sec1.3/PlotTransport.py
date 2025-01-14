import numpy as np
import matplotlib.pyplot as plt

def plotTransport(u, x, t, u0=None):
    """
    Plots the numerical solution of the transport equation and optionally
    compares it with the exact solution derived from the initial condition.

    Parameters
    ----------
    u : ndarray
        Solution array of shape (J+1, N+1), where each column corresponds to time.
    x : ndarray
        Spatial grid points.
    t : ndarray
        Time grid points.
    u0 : function, optional
        Initial condition function u0(x). If provided, the exact solution is plotted for comparison.
    """
    ma = np.max(u)
    mi = np.min(u)
    a = x[0]
    b = x[-1]
    
    # High-resolution spatial grid for exact solution
    xx = np.linspace(a, b, 500)
    
    for n in range(u.shape[1]):  # Iterate over time steps
        plt.clf()
        if u0 is not None:
            exact_solution = u0(xx - t[n])
            plt.plot(x, u[:, n], 'o', label='Numerical Solution', markersize=8)
            plt.plot(xx, exact_solution, '-', label='Exact Solution', linewidth=2)
        else:
            plt.plot(x, u[:, n], 'o', label='Numerical Solution', markersize=8)
        
        plt.xlabel('x')
        plt.ylabel('u')
        plt.ylim([mi, ma])
        plt.title(f'Time t = {t[n]:.2f}')
        plt.legend()
      
