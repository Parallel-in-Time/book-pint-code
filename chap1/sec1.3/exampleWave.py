import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 5000  # Number of time steps
M = 50    # Number of spatial intervals
L = 1     # Length of the domain
dx = L / M
T = 10    # Total simulation time
dt = T / N

x = np.linspace(0, L, M + 1)  # Spatial grid
u = np.zeros((M + 1, N + 1))  # Solution array

# Initial conditions
for j in range(len(x)):
    if x[j] < L / 2:
        u[j, 0] = x[j]  # Triangular initial profile
    else:
        u[j, 0] = L - x[j]
    u[j, 1] = u[j, 0]  # Zero initial velocity (ut0 = 0)

# Time-stepping loop
fig = plt.figure("exampleWave")
for n in range(1, N):
    if not plt.fignum_exists("exampleWave"): break
    
    # Update interior points using finite difference method
    u[1:-1, n + 1] = (
        2 * u[1:-1, n]
        - u[1:-1, n - 1]
        + (dt**2 / dx**2) * (u[0:-2, n] - 2 * u[1:-1, n] + u[2:, n])
    )

    # Plot solution at current time step
    plt.clf()
    plt.plot(x, u[:, n], label=f't={(n-1)*dt:.2f}')
    plt.axis([0, 1, -0.6, 0.6])
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('Wave Equation Simulation')
    plt.grid()
    plt.pause(0.01)  # Pause for animation effect

plt.show()
