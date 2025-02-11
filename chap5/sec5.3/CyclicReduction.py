import numpy as np
import scipy.sparse as sp

def cyclicReduction(A, f):
    """
    Performs a cyclic reduction for a lower bidiagonal system.

    This function performs a cyclic reduction for a lower bidiagonal system Ax=b of even size.

    Parameters
    ----------
    A : sparse matrix
        Coefficient matrix of the system.
    f : array
        Right-hand side vector.

    Returns
    -------
    B : sparse matrix
        Reduced coefficient matrix.
    g : array
        Reduced right-hand side vector.
    """
    n = len(f)
    if n % 2 != 0:
        raise ValueError("The length of f must be even.")

    # Extract diagonals
    d = A.diagonal()
    dm = A.diagonal(-1)

    # Indices for odd and even rows
    idOdd = np.arange(0, n-1, 2)
    idEven = np.arange(1, n, 2)

    # Perform cyclic reduction
    dn = -dm[idEven-1] / d[idOdd] * dm[idEven]
    B = sp.diags([d[idEven], dn], [0, -1], shape=(n//2, n//2))

    # Reduce the right-hand side vector
    g = f[idEven] - dm[idEven-1] / d[idOdd] * f[idOdd]

    return B, g

# Example usage:
# Define a lower bidiagonal matrix A and a vector f
A = sp.diags([-np.ones(5), np.ones(6)], [-1, 0], shape=(6, 6))
f = np.array([1, 2, 3, 4, 5, 6])

# Perform cyclic reduction
B, g = cyclicReduction(A, f)

# Display the results
print("Reduced matrix B:")
print(B.toarray())
print("Reduced vector g:")
print(g)
