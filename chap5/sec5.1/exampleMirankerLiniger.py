import numpy as np
import matplotlib.pyplot as plt

from MirankerLinigerS import mirankerLinigerS
from MirankerLinigerP import mirankerLinigerP

parallel = True

sigma = 10
r = 28
b = 8 / 3
lorenz = lambda t, x: np.array([
    sigma * (x[1] - x[0]),
    r * x[0] - x[1] - x[0] * x[2],
    x[0] * x[1] - b * x[2]
])

# Time parameters
T = 30
N = 30000
dt = T / N
tspan = (0, T)

# Initial conditions
x0 = np.array([20, 5, -5])

# Solve the ODE using MirankerLinigerS method
solver = mirankerLinigerP if parallel else mirankerLinigerS
t, xS = mirankerLinigerP(lorenz, tspan, x0, N)

# Plot the solution
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(xS[0, :], xS[1, :], xS[2, :], '-b')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Lorenz System')
plt.show()
