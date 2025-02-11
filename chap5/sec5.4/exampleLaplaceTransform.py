import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Parameters
J = 100
e = np.ones(J)
dx = 1 / (J + 1)

# Finite difference Laplacian
L = sp.spdiags([e, -2*e, e], [-1, 0, 1], J, J) / dx**2

# Initial temperature
u0 = 20 * e

# Number of quadrature points
N = 100
t = 0.1
ga = np.pi**2 / 2
tau = t / 2
al = tau / (2 * t)

# Quadrature nodes and weights
sigma = 1 / (al * t) * np.log(N / np.arange(1, N+1))
w = np.exp(-ga * t) / (al * t * N) * (np.arange(1, N+1) / N)**(1/al - 1)
w[-1] = w[-1] / 2

# Simpson's rule weights
ws = 2 * w / 3
ws[::2] = ws[::2] * 2

# Function g
def g(t, sigma, uh):
    return 1 / np.pi * np.imag((1 + 1j) * np.exp(-1j * sigma * t) * uh)

# Solution vector
u = np.zeros(J)

# Quadrature loop
for n in range(N):
    s = ga + (1 + 1j) * sigma[n]
    A = L.astype(complex) + s * sp.eye(J)
    uh = -spla.spsolve(A.tocsr(), u0)
    u = u + w[n] * g(t, sigma[n], uh)

# Plot the solution
x = np.linspace(dx, 1-dx, J)
plt.plot(x, u)
plt.xlabel('x')
plt.axis([0, 1, 0, 21])
plt.title('Temperature Distribution')
plt.show()
