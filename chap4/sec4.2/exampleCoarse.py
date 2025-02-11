import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, eye

# Parameters
J = 63                                          # Number of fine grid points
Jc = (J + 1) // 2 - 1                           # Number of coarse grid points

# Prolongation matrix (interpolation)
P = lil_matrix((J, Jc))                         # Initialize sparse matrix
for j in range(Jc):
    P[2*j+1, j] = 1
    P[2*j, j] = 0.5
    P[2*j + 2, j] = 0.5
P = P.tocsr()

# Restriction matrix (transpose of prolongation, scaled by 0.5)
R = 0.5 * P.T

# Discrete Laplacian matrix (fine grid)
dx = 1 / (J + 1)
e = np.ones(J)
L = 1 / dx**2 * csr_matrix(np.diag(-2 * e) + np.diag(e[:-1], k=1) + np.diag(e[:-1], k=-1))

# Coarse matrix by Galerkin projection
Lc = R @ L @ P

# Coarsened time-stepping matrix
dt = 0.01                                       # Time step size
Ac = eye(Lc.shape[0]) - dt * Lc                 # Coarse matrix

# Coarse spatial mesh
xc = np.linspace(0, 1, J + 2)[::2]              # Coarse spatial points
