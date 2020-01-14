import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation


def animate(nframe):
    plt.cla()
    plt.plot(x, f1[nframe][:-1], label=f"Pe = {PE1}")
    plt.plot(x, f2[nframe][:-1], label=f"Pe = {PE2}")
    plt.plot(x, f3[nframe][:-1], label=f"Pe = {PE3}")
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.title(f"Convection-diffusion eq. (t = {nframe})")
    plt.legend(loc='upper right')
    plt.xlabel('x')
    plt.ylabel('f')
    ax.grid(which="both")


# setup condition
PE = 50
D = 1.0e-9
L = 10.0
# u = PE * D / L
dx = 0.05
dt = 625000
# c = u * dt / dx
# d = D * dt / dx / dx
d = 0.2
c = 0.7

t_i = 101
x_i = int(L / dx)
x = np.zeros(x_i + 1)
for i in range(x_i + 1):
    x[i] = dx * i
f = np.zeros((t_i + 1, x_i + 2))

# I.C.
for i in range(int(1.0 / dx) + 1):
    f[0][i] = 0.5 * (1 - np.cos(2 * np.pi * x[i]))

if __name__ == '__main__':
    for i in range(t_i):
        for j in range(1, x_i + 1):
            f[i + 1][j] = f[i][j] - c * 0.5 * (f[i][j + 1] - f[i][j - 1]) + d * (f[i][j + 1] - 2 * f[i][j] + f[i][j - 1])
            if j == x_i - 1:  # Norman condition
                f[i + 1][j + 2] = f[i + 1][j]

    # graph
    fig1, ax1 = plt.subplots()
    p1 = plt.plot(x, f[0][:-1], label="t = 0")
    p2 = plt.plot(x, f[20][:-1], label=f"t = {int(dt * 20)}")
    p3 = plt.plot(x, f[40][:-1], label=f"t = {int(dt * 40)}")
    p4 = plt.plot(x, f[60][:-1], label=f"t = {int(dt * 60)}")
    p5 = plt.plot(x, f[80][:-1], label=f"t = {int(dt * 80)}")
    p6 = plt.plot(x, f[100][:-1], label=f"t = {int(dt * 100)}")
    plt.title("t: 1-order")
    plt.legend(loc='upper right')
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('f')
    # plt.show()

    fig, ax = plt.subplots()
    anim = animation.FuncAnimation(fig, animate, frames=500)
    os.makedirs('../../data/lecture6/images', exist_ok=True)
    anim.save('../../data/convection_diffusion_6_4_1.gif', writer='imagemagick', fps=30)
