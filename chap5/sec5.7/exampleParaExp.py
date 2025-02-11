import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from scipy.linalg import expm

# Parameters
T = 0.4
g = lambda x, t: 10 * np.ones_like(x)
p = 4
N = 10
J = 10
e = np.ones(J-1)
dx = 1 / J
x = np.linspace(dx, 1-dx, J-1)

# Finite difference Laplacian
A = sp.spdiags([e, -2*e, e], [-1, 0, 1], J-1, J-1) / dx**2

# Time steps
dt = T / (p * N)
dT = T / p

# Initial temperature
u0 = np.ones(J-1)
u = np.zeros((J-1, N*p+1))
u[:, 0] = u0

# Compute reference solution
for n in range(N*p):
    u[:, n+1] = spla.spsolve(sp.eye(J-1) - dt * A, u[:, n] + dt * g(x, dt * (n+1)))

# Plot reference solution
fig = plt.figure("exampleParaExp")
ax = fig.add_subplot(111, projection='3d')
X, T_mesh = np.meshgrid(x, np.linspace(0, T, N*p+1))
ax.plot_surface(X, T_mesh, u.T, cmap='viridis',
                rstride=1, cstride=1, shade=False)
ax.set_xlabel('x')
ax.set_ylabel('t')
plt.pause(0.5)

# Compute v solutions
v = np.zeros((p, J-1, N+1))
for j in range(p):
    for n in range(N):
        v[j, :, n+1] = spla.spsolve(sp.eye(J-1) - dt * A,
                                    v[j][:, n] + dt * g(x, dt*(j*N+n+1)))

    # Plot v solutions
    T_mesh_v = np.linspace(j * N * dt, (j + 1) * N * dt, N+1)
    ax.plot_surface(*np.meshgrid(X[0], T_mesh_v), v[j].T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    plt.pause(0.5)

# Compute w solutions
w = np.zeros((p, J-1, p+1))
w[0, :, 0] = u0
for j in range(p):
    for n in range(j, p):
        w[j, :, n+1] = expm(dT * A.toarray()) @ w[j, :, n]

    if j < p-1:
        w[j+1, :, j+1] = v[j, :, -1]

    # Plot w solutions
    T_mesh_w = np.linspace(0, T, p+1)
    ax.plot_surface(*np.meshgrid(X[0], T_mesh_w), w[j].T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    plt.pause(0.5)


# Sum ParaExp solution
UPE = w[0].copy()
ax.plot_surface(*np.meshgrid(X[0], T_mesh_w), UPE.T, cmap='viridis',
                rstride=1, cstride=1, shade=False)
ax.set_xlabel('x')
ax.set_ylabel('t')
plt.pause(0.5)

for j in range(1, p):
    UPE += w[j]

    ax.plot_surface(*np.meshgrid(X[0], T_mesh_w), UPE.T, cmap='viridis',
                    rstride=1, cstride=1, shade=False)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    plt.pause(0.5)

UPE[:, p] += v[p-1][:, N]

ax.plot_surface(*np.meshgrid(X[0], T_mesh_w), UPE.T, cmap='viridis',
                rstride=1, cstride=1, shade=False)
ax.set_xlabel('x')
ax.set_ylabel('t')
plt.pause(0.5)
