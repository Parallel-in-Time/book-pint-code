import math
import numpy as np
import matplotlib.pyplot as plt

from DahlquistBE import dahlquistBE
from Parareal import parareal

# Parameters
la = -1
u0 = 1+0j
T = 1
N = 10
K = 10
MF = 20
MG = 1

# Fine and coarse solvers
F = lambda t0, t1, u0: dahlquistBE(la, [t0, t1], u0, MF)[-1][-1]
G = lambda t0, t1, u0: dahlquistBE(la, [t0, t1], u0, MG)[-1][-1]

# Solve with Parareal and Backward Euler
U = parareal(F, G, T, u0, N, K)
t, u = dahlquistBE(la, [0, T], u0, MF*N)

# Compute errors
err = []
for k in range(K+1):
    err.append(np.max(np.abs(u[::MF] - U[k])))

DT = T / N
R0 = abs(1 / (1 - la*DT))
if R0 < 1:
    R0 = 1

errsup = [err[0]]
errlin = [err[0]]

for k in range(1, K):
    term = abs(np.exp(la * DT) - 1 / (1 - la * DT))**k / math.factorial(k)
    sup_bound = err[0] * term * R0**(N - k - 1) * np.prod(np.arange(N - k + 1, N + 1))
    lin_bound = err[0] * (abs(np.exp(la * DT) - 1 / (1 - la * DT)) / (1 - abs(1 / (1 - la * DT))))**k
    errsup.append(sup_bound)
    errlin.append(lin_bound)

# Plot
plt.semilogy(err, '--', label='Parareal Error')
plt.semilogy(errsup, '-', label='Superlinear Bound')
plt.semilogy(errlin, '-', label='Linear Bound')
plt.xlabel('k')
plt.ylabel('Error')
plt.legend()
plt.grid()
plt.show()
