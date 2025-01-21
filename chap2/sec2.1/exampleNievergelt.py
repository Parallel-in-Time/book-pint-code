import numpy as np
from Nievergelt import nievergelt
from ForwardEuler import forwardEuler

f = lambda t,u: np.cos(t)*u                     # RHS of the ODE problem
u0 = 1; T = 2*np.pi                             # initial value, final time
N = 10                                          # number of subintervals
nSteps = 100                                    # fine steps per subinterval
tPred, uPred = forwardEuler(f, [0, T], u0, N)   # approximate prediction
Mn = 2; width = 0.75                            # algorithm parameters
solver = lambda t0, t1, u0: forwardEuler(f, [t0, t1], u0, nSteps)[-1].ravel()
U, uTraj = nievergelt(solver, T, u0, N, uPred, Mn, width)
