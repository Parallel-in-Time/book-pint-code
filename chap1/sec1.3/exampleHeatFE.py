import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

N, J = 2000, 100
e = np.ones(J-1)
dx = 1/J
A = sps.dia_matrix(([e, -2*e, e], [-1, 0, 1]), shape=(J-1, J-1))/dx**2
T = 1/10
dt = T/N
u = np.zeros((N+1, J-1))
u[0] = 20*e

plt.figure("exampleHeatFE")
for n in range(N):
    if not plt.fignum_exists("exampleHeatFE"): break

    u[n+1] = u[n] + dt*A*u[n]
    plt.cla()
    plt.plot(np.arange(0, 1+dx, dx), [0, *u[n+1], 0])
    plt.title(f"t={(n+1)*dt:1.4f}s")
    plt.grid(True), plt.xlabel("x"), plt.ylabel("y"), plt.ylim(-1, 21)
    plt.pause(0.0001)
plt.show()
