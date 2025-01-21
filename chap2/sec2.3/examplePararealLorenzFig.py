import numpy as np
import matplotlib.pyplot as plt

from ForwardEuler import forwardEuler
from Parareal import parareal

# Lorentz parameters and functions
u0 = [20., 5., -5.]
sigma, rho, beta = 10., 28., 8/3

f = lambda t, u : np.array([
    sigma*(u[1]-u[0]),
    u[0]*(rho-u[2])-u[1],
    u[0]*u[1]-beta*u[2]
    ])

MF = 10; MG = 1                                                     # F and G time steps
F = lambda t0,t1,u0: forwardEuler(f, [t0, t1], u0, MF)[-1][-1]      # fine solver F
G = lambda t0,t1,u0: forwardEuler(f, [t0, t1], u0, MG)[-1][-1]      # coarse solver F
K = 20; u0 = [20, 5, -5]; T = 5; N = 500                            # Parareal parameters
U = parareal(F, G, T, u0, N, K)                                     # solve with Parareal
t, u = forwardEuler(f,[0, T], u0, MF*N);                            # fine solution
TT = np.linspace(0, T, N+1)                                         # coarse time mesh

err = np.zeros(K+1)
for k in range(K):
    err[k] = np.max(np.abs(u[::MF] - U[k]))
plt.semilogy(err)
plt.xlabel("k")
plt.ylabel("error")
plt.tight_layout()
plt.savefig("PararealLorenzError.eps")
