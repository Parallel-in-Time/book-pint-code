import numpy as np
import matplotlib.pyplot as plt

from HeatEquationBE import heatEquationBE
from Parareal import parareal

# Problem setup
f = lambda x, t: x**4 * (1 - x)**4 + 10 * np.sin(8 * t)  # Heat source function
T = 8
N = 16
K = 16
J = 10
u0 = np.zeros(J + 1)

# Coarse solver
MG = 1
gl = np.zeros(MG + 1)
gr = np.zeros(MG + 1)
G = lambda t0, t1, u0: heatEquationBE(f, (t0, t1), (0, 1), u0, gl, gr)[-1][-1]

# Fine solver
MF = 10
gl = np.zeros(MF + 1)
gr = np.zeros(MF + 1)
F = lambda t0, t1, u0: heatEquationBE(f, (t0, t1), (0, 1), u0, gl, gr)[-1][-1]

# Parareal computation
U = parareal(F, G, T, u0, N, K)

# Fine solution
glf = np.zeros(MF * N + 1)
grf = np.zeros(MF * N + 1)
u_fine = heatEquationBE(f, (0, T), (0, 1), u0, glf, grf)[-1]

# Mesh parameters
dt = T / (MF * N)
dT = T / N
dx = 1 / J
t = np.linspace(0, T, MF*N + 1)
TT = np.linspace(0, T, N+1)
x = np.linspace(0, 1, J+1)

# Plotting
fig = plt.figure("examplePararealHeat")
ax = fig.add_subplot(111, projection='3d')

for k in range(K):
    if not plt.fignum_exists("examplePararealHeat"): break

    # Reconstruct fine solution from parareal for plotting
    up = np.zeros_like(u_fine)
    for n in range(N):
        up[n * MF:(n + 1) * MF + 1, :] = heatEquationBE(
            f, (n * dT, (n + 1) * dT), (0, 1), U[k][n, :], gl, gr
        )[-1]
    up[N * MF, :] = U[k, -1, :]

    # Plot parareal approximation
    ax.cla()
    T_grid, X_grid = np.meshgrid(t, x)
    ax.plot_surface(X_grid, T_grid, up.T, cmap='viridis')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title(f'Parareal Approximation (Iteration {k})')
    plt.pause(1)

    # Compute and plot error
    error = np.abs(u_fine - up)
    ax.cla()
    ax.plot_surface(X_grid, T_grid, error.T, cmap='hot')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_title(f'Parareal Error (Iteration {k})')
    plt.pause(1)

plt.show()
