import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Parameters
l = 7
J = 2**l - 1
N = J
e = np.ones(J)
dx = 1 / (J + 1)
x = np.linspace(0, 1, J+2)

# Matrix A
diagonals = np.array([e, -2*e, e]) / dx**2
A = sp.spdiags(diagonals, [-1, 0, 1], J, J)

# Time grid
T = 5
dt = T / N
t = np.linspace(0, T, N+1)

# Initial and boundary conditions
u0 = lambda x: np.zeros_like(x)
gl = lambda t: np.zeros_like(t)
gr = lambda t: np.zeros_like(t)
f = lambda x, t: x**4 * (1-x)**4 + 10 * np.sin(8*t)

# Compute rhs b for reuse
b = np.zeros((J, N))
for n in range(N):
    b[:, n] = dt * f(x[1:-1], t[n+1])
    b[0, n] += dt / dx**2 * gl(t[n+1])
    b[-1, n] += dt / dx**2 * gr(t[n+1])

# Initial solution
u = np.zeros((J+2, N+1))
u[:, 0] = u0(x)
u[0, :] = gl(t)
u[-1, :] = gr(t)

# BE time stepping matrix
G = sp.eye(J) - dt * A

# Compute exact solution
for n in range(N):
    u[1:-1, n+1] = spla.spsolve(G, u[1:-1, n] + b[:, n])

uBE = u.copy()

# Coarse grid setup in space
Jc = (J + 1) // 2 - 1
P = sp.lil_matrix((J, Jc))
for j in range(Jc):
    P[2*j+1, j] = 1
    P[2*j, j] = 0.5
    P[2*j + 2, j] = 0.5

P = P.tocsc()
R = 0.5 * P.T
Ac = R @ A @ P

# Coarse grid setup in time
Nc = (N + 1) // 2 - 1
Pt = sp.lil_matrix((N, Nc))
for j in range(Nc):
    Pt[2*j+1, j] = 1
    Pt[2*j, j] = 0.5
    Pt[2*j + 2, j] = 0.5

Pt = Pt.tocsc()
Rt = 0.5 * Pt.T
Pt[-1, -1] = 1  # No final zero bc in time
Gc = sp.eye(Jc) - 2 * dt * Ac

# Multigrid iteration
nu = 5
al = 0.5
np.random.seed(0)
u[1:-1, 1:] = np.random.rand(J, N)
errcxt = np.zeros(11)
errcxt[0] = np.max(np.abs(uBE - u))

fig = plt.figure("exampleParabolicMultigrid", layout='tight')
ax = fig.add_subplot(111, projection='3d')

for k in range(10):
    for _ in range(nu):  # Block Jacobi steps
        uo = u.copy()
        for n in range(N):
            u[1:-1, n+1] = (1-al) * uo[1:-1, n+1] + al * spla.spsolve(G, uo[1:-1, n] + b[:, n])

    # Plot error after presmoothing
    X, T = np.meshgrid(x, t)

    ax.cla()
    ax.plot_surface(X, T, (uBE - u).T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title('Error after presmoothing')
    plt.pause(0.5)

    # Compute residual
    r = np.zeros((J, N))
    for n in range(N):
        r[:, n] = u[1:-1, n] + b[:, n] - G * u[1:-1, n+1]

    # Restrict residual in space and time
    rc = R * r
    rc = rc @ Rt.T

    # Coarse correction
    uc = np.zeros((Jc+2, Nc+1))
    for n in range(Nc):
        uc[1:-1, n+1] = spla.spsolve(Gc, uc[1:-1, n] + rc[:, n])

    # Extend in space and time
    u[1:-1, 1:] += P @ uc[1:-1, 1:] @ Pt.T

    # Plot error after coarse correction
    ax.cla()
    ax.plot_surface(X, T, (uBE - u).T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title('Error after coarse correction')
    plt.pause(0.5)

    for _ in range(nu):  # Block Jacobi steps
        uo = u.copy()
        for n in range(N):
            u[1:-1, n+1] = (1-al) * uo[1:-1, n+1] + al * spla.spsolve(G, uo[1:-1, n] + b[:, n])

    errcxt[k+1] = np.max(np.abs(uBE - u))

    # Plot error after 2-grid iteration
    ax.cla()
    ax.plot_surface(X, T, (uBE - u).T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title(f'Error after 2-grid iteration k={k}')
    plt.pause(0.5)
