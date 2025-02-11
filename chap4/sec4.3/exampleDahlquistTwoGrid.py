import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, lil_matrix
from scipy.sparse.linalg import inv

# Parameters
l = 6                                   # number of levels
N = 2**l - 1                            # number of gridpoints in time
T = 1
dt = T / N
t = np.linspace(0, T, N + 1)            # time grid
la = -1                                 # Dahlquist parameter

# BE time stepping matrix
e = np.ones(N)
A = diags([-e, (1 - dt * la) * e], [-1, 0], shape=(N, N))

u0 = 0
al = 0.5                                # Jacobi relaxation parameter

# Random initial guess
np.random.seed(0)                       # Set random seed for reproducibility
u = np.random.rand(N)

# Coarse grid size
Nc = 2**(l - 1) - 1

# Prolongation matrix P (sparse)
P = lil_matrix((N, Nc))                 # Initialize sparse matrix
for j in range(Nc):
    P[2 * j, j] = 1
    P[2 * j - 1, j] = 0.5
    P[2 * j + 1, j] = 0.5
P = P.tocsr()                           # Convert to CSR format for efficiency

# Restriction matrix R
R = 0.5 * P.T                           # Transpose of P scaled by 0.5

# Coarse matrix by Galerkin
Ac = R @ A @ P                          # Coarse grid matrix

# Number of presmoothing steps
nu = 4

# Main loop
err = []
for k in range(10):
    err.append(np.max(np.abs(u)))

    # Presmoothing
    for i in range(nu):
        u = u - (al / (1 - dt * la)) * A @ u
        plt.plot(t, np.concatenate(([u0], u)), '-')
        plt.xlabel('t')
        plt.ylabel('error')
        plt.pause(0.1)                   # Pause to visualize the plot

    # Compute coarse correction
    rc = R @ (-A @ u)                    # Restrict residual to coarse grid
    uc = inv(Ac) @ rc                    # Solve on coarse grid
    u = u + P @ uc                       # Prolong and update solution

    # Plot after coarse correction
    plt.plot(t, np.concatenate(([u0], u)), '-r')
    plt.legend(['before coarse', 'after coarse'])
    plt.pause(0.1)                       # Pause to visualize the plot

# Plot error over iterations
plt.figure()
plt.plot(range(1, 11), err, '-o')
plt.xlabel('Iteration')
plt.ylabel('Max Error')
plt.title('Error over Iterations')
plt.show()
