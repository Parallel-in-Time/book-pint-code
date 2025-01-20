import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

# Parameters
l = 4
J = 2**l - 1  # Number of grid points
e = np.ones(J)
h = 1 / (J + 1)
x = np.linspace(0, 1, J + 2)  # Spatial grid including boundaries

# Discrete Laplacian matrix (1D)
A = 1 / h**2 * sp.diags([e, -2 * e, e], [-1, 0, 1], shape=(J, J))

# Damping parameter
w = 1

# Initial random error (consistent for reproducibility)
np.random.seed(0)
u = np.random.rand(J)

# Jacobi iteration

plt.figure("exampleJacobi")
for i in range(20):  # 20 iterations
    if not plt.fignum_exists("exampleJacobi"): break

    plt.cla()
    plt.plot(x, np.concatenate(([0], u, [0])), label=f"Iteration {i}")
    plt.xlabel("x")
    plt.ylabel(f"Error (iteration {i})")
    plt.axis([0, 1, -0.1, 1])
    plt.grid()
    plt.legend()
    plt.pause(0.5)  # Pause for animation
    
    # Jacobi iteration step (error update, f = 0)
    u = u + w * h**2 / 2 * (A @ u)

plt.show()
