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
    iOdd = np.arange(0, n-1, 2)
    iEven = iOdd+1

    # Perform cyclic reduction
    dn = - dm[iOdd[1:]] / d[iOdd[1:]] * dm[iEven[:-1]]
    B = sp.diags([d[iEven], dn], [0, -1], shape=(n//2, n//2))

    # Reduce the right-hand side vector
    g = f[iEven] - dm[iOdd] / d[iOdd] * f[iOdd]

    return B, g
