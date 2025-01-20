import numpy as np
import matplotlib.pyplot as plt

sigma, r, b = 10, 28, 8/3
def f(t, x):
    x1, x2, x3 = x
    return np.array([sigma*(x2-x1), r*x1-x2-x1*x3, x1*x2-b*x3])
T, N = 30, 30000
dt = T / N

x = np.zeros((3, N+1))
x[:, 0] = [20, 5, -5]

ax = plt.figure("exampleLorenz").add_subplot(projection='3d')
plt.grid()

xf = np.sqrt(b * (r - 1))       # Fixed points
yf = np.sqrt(b * (r - 1))
zf = r - 1

for i in range(N):
    if not plt.fignum_exists("exampleLorenz"): break

    x[:, i + 1] = x[:, i] + dt*f(i*dt, x[:, i])     # Forward Euler step
    if (i+1) % 100 == 0:                            # plot only every 100th
        ax.cla()                                    # for animation speed
        ax.plot(x[0, :i + 1], x[1, :i + 1], x[2, :i + 1], '-b', lw=1)
        ax.set_xlabel('x'), ax.set_xlim(-20, 30)
        ax.set_ylabel('y'), ax.set_ylim(-30, 40)
        ax.set_zlabel('z'), ax.set_zlim(-10, 60)
        ax.scatter([0], [0], [0], marker='o', color="red") 
        ax.scatter([xf, -xf], [yf, -yf], [zf, zf], marker='o', color="red")
        ax.set_title(f"t={(i+1)*dt:1.1f}s")
        plt.pause(0.001)

plt.show()
