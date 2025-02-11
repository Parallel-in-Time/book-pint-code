import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags

# Parameters
l = 6                                   # can later become number of levels
N = 2**l - 1                            # number of gridpoints in time
T = 1
dt = T / N
t = np.linspace(0, T, N+1)              # time grid
la = -1                                 # Dahlquist parameter

# BE time stepping matrix
e = np.ones(N)
A = diags([-e, (1 - dt * la) * e], [-1, 0], shape=(N, N))

u0 = 0
al = 0.5                                # Jacobi relaxation parameter

# Random initial guess
np.random.seed(0)                       # Set random seed for reproducibility
u = np.random.rand(N)

# Iteration loop
for i in range(1, 21):
    plt.plot(t, np.concatenate(([u0], u)), '-')
    plt.xlabel('t')
    plt.ylabel(f'error iter = {i-1}')
    plt.show()

    # Damped Jacobi iteration
    u = u - (al / (1 - dt * la)) * A.dot(u)

    plt.pause(0.1)  # Pause to visualize the plot
