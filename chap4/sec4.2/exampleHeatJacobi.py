import numpy as np
from scipy.sparse import spdiags, eye
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

# Parameters
l = 7
J = 2**l - 1                                   # Mesh points in space
e = np.ones(J)
dx = 1 / (J + 1)
x = np.linspace(0, 1, J + 2)                   # Spatial mesh including boundaries
L = 1 / dx**2 * spdiags([e, -2 * e, e], [-1, 0, 1], J, J)  # Discrete Laplacian

T = 5
N = J                                          # Number of time steps
dt = T / N                                     # Time step size
t = np.linspace(0, T, N + 1)                   # Time mesh

# Functions
u0 = lambda x: np.zeros_like(x)                # Initial condition
gl = lambda t: np.zeros_like(t)                # Boundary conditions
gr = lambda t: np.zeros_like(t)
f = lambda x, t: x**4 * (1 - x)**4 + 10 * np.sin(8 * t)  # Source function

# Precompute the right-hand side b
b = np.zeros((J, N))
for n in range(N):
    b[:, n] = dt * f(x[1:-1], t[n + 1])        # Source function
    b[0, n] += dt / dx**2 * gl(t[n + 1])       # Boundary conditions
    b[-1, n] += dt / dx**2 * gr(t[n + 1])

# Initial and boundary conditions
u = np.zeros((J + 2, N + 1))
u[:, 0] = u0(x)                                # Initial condition
u[0, :] = gl(t)                                # Left boundary condition
u[-1, :] = gr(t)                               # Right boundary condition

# Time-stepping matrix
A = eye(J) - dt * L                            # Backward Euler matrix

# Compute exact Backward Euler solution
for n in range(N):
    u[1:-1, n + 1] = spsolve(A, u[1:-1, n] + b[:, n])

# Store the exact solution
uBE = u.copy()

# Jacobi iteration
D = spdiags(np.diag(A.toarray()), 0, J, J).tocsr()      # Diagonal part of A
NU = [5000, 1000, 100, 5]                               # Number of Jacobi steps

fig = plt.figure("exampleHeatJacobi")
ax = fig.add_subplot(111, projection='3d')

for nu in NU:
    if not plt.fignum_exists("exampleHeatJacobi"): break

    u = uBE.copy()                             # Reset to exact solution
    for n in range(N):
        v = u[1:-1, n]                         # Initial guess from the previous time step
        print(f"{n+1}/{N}")
        for _ in range(nu):
            v = v + spsolve(D, (u[1:-1, n] + b[:, n] - A @ v))  # Jacobi iteration
        u[1:-1, n + 1] = v                     # Update the solution

    # Plot the results
    try:
        ax.collections[-1].colorbar.remove()
    except: pass
    ax.cla()
    X, T = np.meshgrid(x, t)
    c = ax.plot_surface(X, T, u.T, cmap='viridis', rstride=1, cstride=1)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title(f'Approximation after nu = {nu} Jacobi steps')
    ax.set_zlim(-1, 1)
    plt.colorbar(c, ax=ax)
    plt.pause(0.5)

plt.show()
