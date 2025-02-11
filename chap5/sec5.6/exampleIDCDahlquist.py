import numpy as np
import matplotlib.pyplot as plt
from RIDC import RIDC

# Parameters
lam = 1j - 0.02
T = 6 * np.pi
M = 5
N = 10
u0 = 1

# Solve using IDC with K=0
t, uIDC = RIDC(lam, [0, T], u0, N, M, 0)

# Exact solution
uExact = np.exp(lam * t)

# Plot the results
plt.plot(t, np.real(uExact), 'o-', label='Exact Solution')
plt.plot(t, np.real(uIDC), label='IDC K=0')

# Solve and plot for K=1 to K=3
for K in range(1, 4):
    t, uIDC = RIDC(lam, [0, T], u0, N, M, K)
    plt.plot(t, np.real(uIDC), label=f'IDC K={K}')

plt.xlabel('t')
plt.ylabel('Real(u)')
plt.title('IDC Solutions vs Exact Solution')
plt.legend()
plt.show()
