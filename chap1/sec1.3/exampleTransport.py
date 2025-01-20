import numpy as np

from TransportFEUpwind import transportFEUpwind
from PlotTransport import plotTransport

f = lambda x, t: np.zeros_like(x)               # Zero right-hand side
u0 = lambda x: np.exp(-120 * (x - 0.3) ** 2)    # Initial condition
g = lambda t: np.zeros_like(t)                  # Boundary condition

a, b = 0, 1     # Spatial domain
J = 40          # Number of spatial intervals
T = 0.5         # Final time
N = 20          # Number of time steps

u, x, t = transportFEUpwind(f, u0, g, a, b, J, T, N)
plotTransport(u, x, t, u0)
