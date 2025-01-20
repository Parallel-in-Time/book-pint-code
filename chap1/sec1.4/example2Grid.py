import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Parameters
l = 4
J = 2**l - 1  # Fine grid size
e = np.ones(J)
h = 1 / (J + 1)
x = np.linspace(0, 1, J + 2)

# Construct the 1D Laplacian matrix on the fine grid
A = 1 / h**2 * sp.diags([e, -2 * e, e], [-1, 0, 1], shape=(J, J))

# Coarse grid parameters
Jc = 2**(l - 1) - 1  # Coarse grid size
P = sp.lil_matrix((J, Jc))  # Prolongation matrix (interpolation)
for j in range(1, Jc + 1):
    P[2 * j - 1, j - 1] = 1
    if 2 * j - 2 >= 0:
        P[2 * j - 2, j - 1] = 0.5
    if 2 * j < J:
        P[2 * j, j - 1] = 0.5
P = P.tocsr()

R = 0.5 * P.T  # Restriction matrix
Ac = R @ A @ P  # Galerkin coarse-grid matrix

# Multigrid parameters
nu = 2  # Number of smoothing steps
np.random.seed(0)
u = np.random.rand(J)  # Random initial guess
w = 2 / 3  # Jacobi damping parameter

# Multigrid V-cycle
plt.figure("example2Grid")
for n in range(3):  # Number of multigrid iterations
    if not plt.fignum_exists("example2Grid"): break

    for i in range(nu):  # Jacobi smoothing steps
        u = u + w * h**2 / 2 * A @ u  # Jacobi damping step

    # Plot before coarse correction
    plt.cla()
    plt.plot(x, np.concatenate(([0], u, [0])), label="Before coarse correction")
    plt.xlabel("x")
    plt.ylabel(f"Error (iteration {n + 1})")
    plt.axis([0, 1, -0.1, 1])
    
    # Compute residual and coarse correction
    rc = R @ (-A @ u)  # Residual (f = 0)
    uc = spla.spsolve(Ac, rc)  # Solve coarse problem
    u = u + P @ uc  # Correct fine-grid solution

    # Plot after coarse correction
    plt.plot(x, np.concatenate(([0], u, [0])), "r", label="After coarse correction")
    plt.legend()
    plt.grid()
    plt.pause(1.5)

plt.show()
