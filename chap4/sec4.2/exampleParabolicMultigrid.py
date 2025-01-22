import numpy as np
from scipy.sparse import spdiags, eye, lil_matrix
import matplotlib.pyplot as plt

# Parameters
l = 7
J = 2**l - 1                                   # Number of spatial mesh points
dx = 1 / (J + 1)
x = np.linspace(0, 1, J + 2)                   # Spatial mesh
e = np.ones(J)
L = 1 / dx**2 * spdiags([e, -2*e, e], [-1, 0, 1], J, J)  # Discrete Laplacian

T = 5
N = J                                          # Number of time steps
dt = T / N
t = np.linspace(0, T, N + 1)                   # Time mesh

# Initial and boundary conditions
u0 = lambda x: np.zeros_like(x)                # Initial condition
gl = lambda t: np.zeros_like(t)                # Left boundary
gr = lambda t: np.zeros_like(t)                # Right boundary

# Source function
f = lambda x, t: x**4 * (1 - x)**4 + 10 * np.sin(8 * t)

# Precompute the RHS (b)
b = np.zeros((J, N))
for n in range(N):
    b[:, n] = dt * f(x[1:-1], t[n + 1])        # Source function
    b[0, n] += dt / dx**2 * gl(t[n + 1])       # Boundary conditions
    b[-1, n] += dt / dx**2 * gr(t[n + 1])

# Initialize solution array
u = np.zeros((J + 2, N + 1))
u[:, 0] = u0(x)                                # Initial condition
u[0, :] = gl(t)                                # Boundary conditions
u[-1, :] = gr(t)

# Time-stepping matrix
A = eye(J) - dt * L
D_inv = 1 / A.diagonal()                       # Inverse of the diagonal of A

# Exact Backward Euler solution
for n in range(N):
    u[1:-1, n + 1] = np.linalg.solve(A.toarray(), u[1:-1, n] + b[:, n])
uBE = u.copy()                                 # Store the exact BE solution

# Multigrid Parameters
Jc = (J + 1) // 2 - 1                          # Coarse grid points
P = lil_matrix((J, Jc))                        # Prolongation matrix
for j in range(Jc):
    if 2 * j - 1 >= 0:
        P[2 * j - 1, j] = 0.5
    P[2 * j, j] = 1
    if 2 * j + 1 < J:
        P[2 * j + 1, j] = 0.5

R = 0.5 * P.T                                  # Restriction matrix
Lc = R @ L @ P                                 # Coarse matrix
Ac = eye(Jc) - dt * Lc                         # Coarse time-stepping matrix

# Multigrid Solver
nu = 5
alpha = 0.5
K = 10                                         # Number of iterations
u = uBE.copy()                                 # Random initial guess
u[1:-1, 1:] = np.random.rand(J, N)             # Random interior guess

fig = plt.figure("exampleParabolicMultigrid")
ax = fig.add_subplot(projection='3d')

for k in range(K):
    if not plt.fignum_exists("exampleParabolicMultigrid"): break

    # Presmoothing
    for n in range(N):
        v = u[1:-1, n + 1]
        for j in range(nu):
            v += alpha * (D_inv * (u[1:-1, n] + b[:, n] - A @ v))  # Correct Jacobi iteration
        u[1:-1, n + 1] = v

    # Plot error after presmoothing
    try:
        ax.collections[-1].colorbar.remove()
    except: pass
    ax.cla()
    X, T = np.meshgrid(x, t)
    c = ax.plot_surface(X, T, uBE.T - u.T, cmap='viridis', rstride=1, cstride=1)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title(f'Error after presmoothing, iteration k={k+1}')
    ax.view_init(15, -140)
    plt.colorbar(c, ax=ax)
    plt.pause(1)
    if not plt.fignum_exists("exampleParabolicMultigrid"): break

    # Compute residual
    r = np.zeros_like(b)
    for n in range(N):
        r[:, n] = u[1:-1, n] + b[:, n] - A @ u[1:-1, n + 1]

    # Coarse grid correction
    rc = R @ r
    uc = np.zeros((Jc + 2, N + 1))             # Zero initial guess
    for n in range(N):
        uc[1:-1, n + 1] = np.linalg.solve(Ac.toarray(), uc[1:-1, n] + rc[:, n])
    u[1:-1, :] += P @ uc[1:-1, :]              # Add coarse correction

    # Plot error after coarse correction
    try:
        ax.collections[-1].colorbar.remove()
    except: pass
    ax.cla()
    X, T = np.meshgrid(x, t)
    c = ax.plot_surface(X, T, uBE.T - u.T, cmap='viridis', rstride=1, cstride=1)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title(f'Error after correction, iteration k={k+1}')
    ax.view_init(15, -140)
    plt.colorbar(c, ax=ax)
    plt.pause(1)
    if not plt.fignum_exists("exampleParabolicMultigrid"): break

    aaaaa

    # Postsmoothing
    for n in range(N):
        v = u[1:-1, n + 1]
        for j in range(nu):
            v += alpha * (D_inv * (u[1:-1, n] + b[:, n] - A @ v))  # Correct Jacobi iteration
        u[1:-1, n + 1] = v

    # Plot error after postsmoothing
    try:
        ax.collections[-1].colorbar.remove()
    except: pass
    ax.cla()
    X, T = np.meshgrid(x, t)
    c = ax.plot_surface(X, T, uBE.T - u.T, cmap='viridis', rstride=1, cstride=1)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title(f'Error after postsmoothing, iteration k={k+1}')
    ax.view_init(15, -140)
    plt.colorbar(c, ax=ax)
    plt.pause(1)
    if not plt.fignum_exists("exampleParabolicMultigrid"): break
