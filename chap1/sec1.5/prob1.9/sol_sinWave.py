#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 15:47:50 2018

@author: lunet
"""
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt
from solver import waveEquation

L = np.pi
c = 1
J = 15
tSpan = 16*np.pi

x = np.linspace(0, L, num=J+2)

u0 = np.sin(2*x)
ut0 = 0

nStep = 300
t1, u1 = waveEquation(u0, ut0, c, L, tSpan, nStep)

nStep = 600
t2, u2 = waveEquation(u0, ut0, c, L, tSpan, nStep)

nStep = 1200
t3, u3 = waveEquation(u0, ut0, c, L, tSpan, nStep)

fig = plt.figure('3D')
ax = p3.Axes3D(fig)
tMesh, xMesh = np.meshgrid(t1, x)
ax.plot_surface(xMesh[:, :40], tMesh[:, :40], u1[:, :40],
                rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u(x,t)')

plt.figure('u(pi/4,T)')
plt.plot(t1, u1[4, :], '^-', label='nStep=300', markevery=0.1)
plt.plot(t2, u2[4, :], 's-', label='nStep=600', markevery=0.1)
plt.plot(t3, u3[4, :], 'o-', label='nStep=1200', markevery=0.1)
uTh = np.sin(2*np.pi/4)*np.cos(2*c*t2)
plt.plot(t2, uTh, '--', label='Analytical')
plt.grid(True)
plt.legend(loc='center right')
plt.xlim(14*np.pi, 53)
plt.xlabel('$t$')
plt.ylabel('$u(\\pi/4,t)$')
plt.savefig('waveSinSol.pdf', bbox_inches='tight')
plt.show()
