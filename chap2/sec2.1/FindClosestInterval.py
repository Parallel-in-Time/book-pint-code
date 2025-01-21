def findClosestInterval(x,y):
    """
    Finds the closest interval in a vector to a number.

    Parameters
    ----------
    x : float
        The value for which to find the closest interval.
    y : list or numpy.ndarray
        A sorted vector defining the intervals.

    Returns
    -------
    int
        The lower index `m` of the interval in `y` that contains `x`.
        If `x` is outside the range of `y`, returns the closest interval.
    """
    n = len(y)
    i = 0
    while i < n and x > y[i]:
        i += 1
    if i == 0:
        return 0
    elif i == n:
        return n-1
    else:
        return i-1