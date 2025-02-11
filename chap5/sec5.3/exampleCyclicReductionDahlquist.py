import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

from CyclicReduction import cyclicReduction
from CyclicBackSubstitution import cyclicBackSubstitution

# Parameters
m = 4
n = 2**m
x0 = 1
lambda_ = -1
dt = 1 / n

# Solve Dahlquist equation
d = np.ones(n)
dm = -(1 + lambda_ * dt) * d
A = [sp.diags([dm, d], [-1, 0], shape=(n, n))]
f = [np.zeros(n)]
f[0][0] = x0 * (1 + lambda_ * dt)

# Cyclic reduction
for i in range(m):
    A.append(cyclicReduction(A[i], f[i])[0])
    f.append(cyclicReduction(A[i], f[i])[1])

# Solve smallest system
x = [None] * (m + 1)
x[m] = np.linalg.solve(A[m].toarray(), f[m])

# Cyclic back-substitution
for i in range(m-1, -1, -1):
    x[i] = cyclicBackSubstitution(A[i], f[i], x[i+1])

# Compare with exact solution
t = np.linspace(0, 1, n+1)
exact_solution = np.exp(lambda_ * t)

# Plot the results
plt.plot(t, np.concatenate(([x0], x[0])), '--', label='Cyclic Reduction')
plt.plot(t, exact_solution, '-', label='Exact Solution')
plt.xlabel('t')
plt.legend()
plt.title('Cyclic Reduction vs Exact Solution')
plt.show()