import os

import numpy as np
from matplotlib import pyplot as plt


# setup condition
L = 10.0
dx = 0.05
dt = 625000
c = 0.5

t_i = 1001
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
            f[i + 1][j] = f[i][j] - c * (f[i][j] - f[i][j - 1])
            if j == x_i - 1:  # Norman condition
                f[i + 1][j + 2] = f[i + 1][j]

    # graph
    fig, ax = plt.subplots()
    p1 = plt.plot(x, f[0][:-1], label="t = 0")
    p2 = plt.plot(x, f[40][:-1], label=f"t = {int(dt * 40)}")
    p3 = plt.plot(x, f[80][:-1], label=f"t = {int(dt * 80)}")
    p4 = plt.plot(x, f[120][:-1], label=f"t = {int(dt * 120)}")
    p5 = plt.plot(x, f[160][:-1], label=f"t = {int(dt * 160)}")
    p6 = plt.plot(x, f[200][:-1], label=f"t = {int(dt * 200)}")
    plt.title(f"c = {c}")
    plt.legend(loc='upper right')
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('f')
    # plt.show()
    os.makedirs('../../data/lecture6/images', exist_ok=True)
    fig.savefig('../../data/lecture6/images/6_2.png')
