import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Parameters
N = 10
e = np.ones(N)
T = 1
dt = T / N
t = np.linspace(0, T, N+1)
la = -1
u0 = 1

# Time stepping matrix
A = sp.spdiags([-e/dt, e/dt-la], [-1, 0], N, N)

# Right-hand side vector
f = np.zeros(N)
f[0] = u0 / dt

# Exact solution
ue = spla.spsolve(A.tocsr(), f)

# Alpha-cyclic matrix
al = 0.1
At = A.tolil()
At[0, -1] = al * At[1, 0]

# Diagonal scaling
D = sp.spdiags(al**(np.arange(N) / N), 0, N, N)

# Eigenvalues
L = sp.spdiags(np.conj(np.fft.fft([0] + [1] + [0]*(N-2))), 0, N, N)

# Initial guess
u = np.zeros(N)
up = u.copy()

# Number of iterations
K = 5

# Iterative solver
err = np.zeros(K+1)
errp = np.zeros(K+1)

for k in range(K+1):
    plt.gca().cla()
    plt.plot(t, np.concatenate(([u0], ue.real)), '-', label='Solution')
    plt.plot(t, np.concatenate(([u0], u.real)), 'o', label='Direct Solve')
    plt.plot(t, np.concatenate(([u0], up.real)), '*', label='ParaDiag Solve')
    plt.xlabel('t')
    plt.legend()
    plt.pause(1)

    # Compute errors
    err[k] = np.linalg.norm(ue - u)
    errp[k] = np.linalg.norm(ue - up)

    # Update solutions
    u = u + spla.spsolve(At.tocsr(), f - A @ u)
    up = up + dt * al**(-1/N) * np.linalg.inv(D.toarray()) @ np.fft.fft(
        np.linalg.solve(al**(-1/N) * (1 - la * dt) * np.eye(N) - L.toarray(),
                        np.fft.ifft(D.toarray() @ (f - A @ up)))
    )
