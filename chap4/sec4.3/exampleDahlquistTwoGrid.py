import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

# Parameters
l = 6
N = 2**l - 1
T = 1
dt = T / N
t = np.linspace(0, T, N+1)
la = -1
e = np.ones(N)

# BE time stepping matrix
diagonals = np.array([-e, (1 - dt * la) * e])
A = sp.spdiags(diagonals, [-1, 0], N, N)

# Initial conditions
u0 = 0
al = 0.5
np.random.seed(0)  # For reproducibility
u = np.random.rand(N)

# Coarse grid setup
Nc = 2**(l-1) - 1
P = sp.lil_matrix((N, Nc))
for j in range(Nc):
    P[2*j+1, j] = 1
    P[2*j, j] = 0.5
    P[2*j + 2, j] = 0.5

P = P.tocsc()
R = 0.5 * P.T
Ac = R * A * P

# Multigrid iteration
nu = 4
err = np.zeros(10)

fig = plt.figure("exampleDahlquistTwoGrid")

for k in range(10):
    if not plt.fignum_exists("exampleDahlquistTwoGrid"): break

    err[k] = max(abs(u))
    for i in range(nu):  # Presmoothing
        u = u - al / (1 - dt * la) * A @ u

    plt.gca().cla()
    plt.plot(t, np.concatenate(([u0], u)), '-')
    plt.xlabel('t')
    plt.ylabel('error')

    # Coarse grid correction
    rc = R @ (-A @ u)
    u = u + P @ sp.linalg.spsolve(Ac, rc)

    # Plot after coarse correction
    plt.plot(t, np.concatenate(([u0], u)), '-r')
    plt.legend(['before coarse', 'after coarse'])
    plt.pause(0.5)
