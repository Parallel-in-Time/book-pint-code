import numpy as np
import scipy.sparse as sp

def cyclicBackSubstitution(A, f, xr):
    """
    Performs back-substitution after cyclic reduction for a bidiagonal system.

    This function performs the back-substitution in the bidiagonal system Ax=b
    when the even unknowns are already computed and given in xr.

    Parameters
    ----------
    A : sparse matrix
        Coefficient matrix of the system.
    f : array
        Right-hand side vector.
    xr : array
        Known solution values for even indices.

    Returns
    -------
    x : array
        Complete solution vector.
    """
    n = len(f)
    if n % 2 != 0:
        raise ValueError("The length of f must be even.")

    iOdd = np.arange(0, n-1, 2)
    x = np.zeros(n)
    x[iOdd+1] = xr  # Insert known solution values

    # Extract needed diagonal entries
    d = A.diagonal()
    dm = A.diagonal(-1)

    # Compute odd solution values
    x[iOdd] = (f[iOdd+1] - d[iOdd+1] * xr) / dm[iOdd]

    return x

# Example usage:
# Define a lower bidiagonal matrix A, a vector f, and known solution values xr
A = sp.diags([-np.ones(5), np.ones(6)], [-1, 0], shape=(6, 6))
f = np.array([1, 2, 3, 4, 5, 6])
xr = np.array([2, 4, 6])

# Perform cyclic back-substitution
x = cyclicBackSubstitution(A, f, xr)

# Display the result
print("Solution vector x:")
print(x)
