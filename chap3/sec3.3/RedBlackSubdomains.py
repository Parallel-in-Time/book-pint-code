import numpy as np
from scipy.sparse import eye

def redBlackSubdomains(n, m, nx):
    """
    Computes a red-black space-time decomposition for a 1D wave equation.

    Parameters
    ----------
    n : int
        Number of spatial steps.
    m : int
        Number of time steps.
    nx : int
        Number of red spatial subdomains.

    Returns
    -------
    Rr : list of lists of sparse matrices
        Red restriction matrices.
    Rb : list of lists of sparse matrices
        Black restriction matrices.
    mtr : int
        Number of red subdomains in time.
    mtb : int
        Number of black subdomains in time.
    """
    idr = np.round(n/nx*np.arange(nx+1)).astype(int)            # Red interface locations
    sxr = [np.arange(idr[i], idr[i+1]) for i in range(nx)]      # Red subdomain space indices

    idb = np.round(n/nx*(np.arange(0.5, nx))).astype(int)       # Black interface locations
    sxb = [np.arange(idb[i], idb[i+1]) for i in range(nx-1)]    # Black subdomain space indices

    mx = max(np.diff(idr).max(), np.diff(idb).max())            # Maximum subdomain width
    mst = min(mx//2 - 1, m)                                     # Maximum tent height (discrete CFL)

    # Time indices for red subdomains
    str_ = {1: np.arange(mst)}
    mtr = 1
    id = mst
    while id < m:
        mtr += 1
        idn = min(id + 2 * mst, m)
        str_[mtr] = np.arange(id, idn)
        id = idn

    # Time indices for black subdomains
    stb = {1: np.arange(min(2 * mst, m))}
    mtb = 1
    id = min(2 * mst, m)
    while id < m:
        mtb += 1
        idn = min(id + 2 * mst, m)
        stb[mtb] = np.arange(id, idn)
        id = idn

    Id = eye(n*m, format='csr')  # To extract the R matrices
    G = np.arange(1, n*m + 1).reshape(m, n)  # Form unknown enumeration

    # Extract red subdomains
    Rr = {}
    for i in range(nx):
        Rr[i] = {}
        for j in range(1, mtr + 1):
            id = G[str_[j], :][:, sxr[i]]
            Rr[i][j] = Id[id.flatten() - 1, :]

    # Extract black subdomains
    Rb = {}
    for i in range(nx - 1):
        Rb[i] = {}
        for j in range(1, mtb + 1):
            id = G[stb[j], :][:, sxb[i]]
            Rb[i][j] = Id[id.flatten() - 1, :]

    return Rr, Rb, mtr, mtb
