import numpy as np
import matplotlib.pyplot as plt
from Nievergelt import nievergelt
from ForwardEuler import forwardEuler

# RHS of the ODE problem
f = lambda t, u: np.cos(t) * u
u0 = 1  # Initial solution
T = 2 * np.pi  # Final time
N = 10  # Number of subintervals
nSteps = 100  # Fine steps per subinterval

# Approximate prediction using Forward Euler
tPred, uPred = forwardEuler(f, [0, T], u0, N)

# Nievergelt's method parameters
Mn = 2
width = 0.75

# Solver function for Nievergelt's method
solver = lambda t0, t1, u0: forwardEuler(f, [t0, t1], u0, nSteps)[-1].ravel()

# Nievergelt's method
U, uTraj = nievergelt(solver, T, u0, N, uPred, Mn, width)

# Fine solution for comparison
tFine, uFine = forwardEuler(f, [0, T], u0, N * nSteps)

# Plotting
plt.plot(tFine, uFine, label='Fine', zorder=3)
plt.xlabel('Time')
plt.ylabel('Solution')
plt.plot(tPred, uPred, 'o', label='Prediction', zorder=2)
plt.plot(tPred, U, '^', label='Nievergelt', zorder=2)

# Plot trajectories
TT = np.linspace(0, T, N + 1)
x = uTraj[0, 0, :]
plt.plot(np.linspace(TT[0], TT[1], nSteps + 1), x, 'k--', label='Trajectories', zorder=1)

for n in range(1, N):
    for m in range(Mn):
        x = uTraj[n, m, :]
        plt.plot(np.linspace(TT[n], TT[n + 1], nSteps + 1), x, 'k--', zorder=1, alpha=0.7)

# Add legend and finalize plot
plt.legend(loc='upper right')
plt.gca().tick_params(labelsize=12)
plt.savefig('NievergeltExampleLinear.eps', format='eps')
plt.show()
