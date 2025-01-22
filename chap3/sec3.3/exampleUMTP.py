import numpy as np
from scipy.sparse import spdiags, kron, eye
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

from ReadBlackSubdomains import redBlackSubdomains

# Parameters
c = 1                                         # Wave speed
n = 200                                       # Number of spatial points
X = 1                                         # Space domain (0, X)
dx = X / (n + 1)                              # Spatial mesh size
ex = np.ones(n)
DDx = 1 / dx**2 * spdiags([ex, -2 * ex, ex], [-1, 0, 1], n, n)  # Second derivative in space

T = 1                                         # Time domain (0, T)
dtOnCFL = dx / c                              # dt on CFL
m = int(np.ceil(T / dtOnCFL))                 # Number of time steps
dt = T / m                                    # Time step size
et = np.ones(m)
DDt = 1 / (c**2 * dt**2) * spdiags([et, -2 * et, et], [-2, -1, 0], m, m)  # Second derivative in time

# All-at-once matrix
A = kron(DDt, eye(n)) - kron(spdiags(et, -1, m, m), DDx)

# Time and space mesh
t = np.arange(dt, T + dt, dt)
x = np.arange(dx, X, dx)

# Initial conditions
u0 = np.exp(-200 * (x - 0.5)**2)
u0t = np.zeros_like(x)
f = np.zeros(A.shape[0])

# Add initial conditions to the right-hand side
f[:n] = u0t / dt + u0 / dt**2 + 0.5 * DDx @ u0
f[n:2 * n] = -u0 / dt**2

# Exact solution
ue = spsolve(A, f)
Ue = ue.reshape((n, m))

# Plot the exact solution
fig = plt.figure("exampleUMTP_exactSolution")
ax = fig.add_subplot(projection='3d')
X, T = np.meshgrid(t, x)
ax.plot_surface(T, X, Ue, cmap='viridis')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_title('Exact Solution')
plt.colorbar(ax.plot_surface(T, X, Ue, cmap='viridis'), ax=ax)
plt.show()

# Space-time decomposition
nx = 5                                        # Number of red subdomains in space
Rr, Rb, mtr, mtb = redBlackSubdomains(n, m, nx)

# Random initial guess
u = np.random.rand(m * n) - 0.5

fig = plt.figure("exampleUMTP")
ax = fig.add_subplot(111, projection='3d')

# Unmapped tent pitching
for jr in range(1, mtr + 1):
    if not plt.fignum_exists("exampleUMTP"): break

    for ir in range(nx):  # Red subdomain solves
        Rr_ir_jr = Rr[ir][jr]
        u += Rr_ir_jr.T @ spsolve(Rr_ir_jr @ A @ Rr_ir_jr.T, Rr_ir_jr @ (f - A @ u))

    # Plot the error
    U = (ue - u).reshape((n, m))


    try:
        ax.collections[-1].colorbar.remove()
    except: pass
    ax.clear()
    ax.plot_surface(T, X, U, cmap='viridis')
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    ax.set_title('Error After Red Solve')
    plt.colorbar(ax.plot_surface(T, X, U, cmap='viridis'), ax=ax)
    plt.pause(1)

    if jr <= mtb:  # Black subdomain solves
        for ib in range(nx - 1):
            Rb_ib_jr = Rb[ib][jr]
            u += Rb_ib_jr.T @ spsolve(Rb_ib_jr @ A @ Rb_ib_jr.T, Rb_ib_jr @ (f - A @ u))

        # Plot the error
        U = (ue - u).reshape((n, m))
        try:
            ax.collections[-1].colorbar.remove()
        except: pass
        ax.clear()
        ax.plot_surface(T, X, U, cmap='viridis')
        ax.set_xlabel('t')
        ax.set_ylabel('x')
        ax.set_title('Error After Black Solve')
        plt.colorbar(ax.plot_surface(T, X, U, cmap='viridis'), ax=ax)
        plt.pause(1)
