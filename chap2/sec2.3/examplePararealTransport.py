import numpy as np
import matplotlib.pyplot as plt

from TransportBE import transportBE
from Parareal import parareal


# Problem setup
f = lambda x, t: 0
a = 1
T = 4
N = 16
K = 16
J = 20

dx = 1 / J
x = np.linspace(0, 1, J + 1)  # Spatial mesh
u0 = np.sin(2 * np.pi * x)    # Initial condition

# Coarse solver G
MG = 1
G = lambda t0, t1, u0: transportBE(f, a, [t0, t1], [0, 1], u0, MG)[-1][-1]

# Fine solver F
MF = 20
F = lambda t0, t1, u0: transportBE(f, a, [t0, t1], [0, 1], u0, MF)[-1][-1]

# Apply Parareal algorithm
U = parareal(F, G, T, u0, N, K)

# Compute fine solution
u = transportBE(f, a, [0, T], [0, 1], u0, N * MF)

# Time parameters
dt = T / (MF * N)
dT = T / N
t = np.linspace(0, T, MF * N + 1)
TT = np.linspace(0, T, N + 1)


# Visualization and error computation
fig = plt.figure("examplePararealTransport")
ax = fig.add_subplot(111, projection='3d')

err = np.zeros(K)
for k in range(K):
    if not plt.fignum_exists("examplePararealTransport"): break

    up = np.zeros_like(u)
    for n in range(N):  # Reconstruct fine solution from parareal approximation
        up[n * MF:(n + 1) * MF + 1, :] = transportBE(
            f, a, [n * dT, (n + 1) * dT], [0, 1], U[k][n, :], MF
        )
    up[N * MF, :] = U[k][-1, :]  # Final value

    # Plot parareal approximation
    ax.cla()
    ax = fig.add_subplot(111, projection='3d')
    T_grid, X_grid = np.meshgrid(t, x, indexing='ij')
    ax.plot_surface(X_grid, T_grid, up, cmap='viridis', rstride=1, cstride=1)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    plt.title(f'Parareal Approximation (k = {k + 1})')
    plt.pause(1)

    # Plot parareal error
    ax.cla()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X_grid, T_grid, u - up, cmap='viridis', rstride=1, cstride=1)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    plt.title(f'Parareal Error (k = {k + 1})')
    plt.pause(1)

    # Compute maximum error
    err[k] = np.max(np.abs(u - up))

plt.show()
