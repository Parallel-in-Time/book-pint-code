import numpy as np

from ForwardEuler import forwardEuler
from MultipleShooting import multipleShooting

sigma=10; r=28; b=8/3                   # Lorenz rhs and Jacobian
f = lambda t, u: np.array([
    sigma*(u[1]-u[0]), r*u[0]-u[1]-u[0]*u[2], u[0]*u[1]-b*u[2]])
jac = lambda t, u: np.array([
    [-sigma, sigma, 0    ],
    [r-u[2], -1   , -u[0]],
    [u[1]  , u[0] , -b   ]])
M = 10

def rhsFull(t, u):          # RHS of the coupled system (f and Jacobian)
    u, V = u[:3], u[3:].reshape((3, 3))
    uEval = f(t, u)
    jacEval = jac(t, u).dot(V)
    return np.ravel([uEval, *jacEval])

def prop(t0, t1, u0):               # propagator for the coupled system (f and Jacobian)
    u0 = np.ravel([u0, *np.eye(3)]) # initial Jacobian is identity
    t, u = forwardEuler(rhsFull, [t0, t1], u0, M)
    return u[-1, :3], u[-1, 3:].reshape((3, 3))

T=1; u0=[20,5,-5]; K=9; N=500                             # multiple shooting parameters
_, uPred = forwardEuler(f, [0, T], u0, N)                   # compute initial guess using N
                                                          # steps of Forward Euler
times, U = multipleShooting(prop, 0, T, u0, N, K, uPred)  # solve with multiple shooting
