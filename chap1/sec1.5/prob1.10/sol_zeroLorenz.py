# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 08:44:22 2018

@author: telu
"""
import numpy as np
from solver import newton
import matplotlib.pyplot as plt

# Set Lorenz system parameters and define functions
sigma, rho, beta = 10, 28, 8/3


def f(u):
    x, y, z = u
    return sigma*(y-x), x*(rho-z)-y, x*y-beta*z


def fA(u):
    x, y, z = u
    mJac = np.array([[-sigma, sigma, 0],
                     [rho-z, -1, -x],
                     [y, x, -beta]])
    return mJac


# %% Apply Newton methods
x0Goal = np.array([(beta*(rho-1))**0.5, (beta*(rho-1))**0.5, rho-1])
x0 = x0Goal - 5
eps, x = newton(f, fA, x0, 10)

# %% Plot error
plt.semilogy(eps, '-o')
plt.xlabel('$N_{iter}$')
plt.ylabel('Error')
plt.grid(True)
plt.savefig('newtonLorentz.pdf', bbox_inches='tight')
print(x0Goal-x)
