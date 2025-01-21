import numpy as np
from ForwardEuler import forwardEuler
N = 1000                                        # number of time steps
f = lambda t,u: np.cos(t)*u                     # RHS function
u0 = 1                                          # initial value
t, u = forwardEuler(f, [0, 2*np.pi], u0, N);    # numerical solution