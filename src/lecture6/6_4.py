import numpy as np
from matplotlib import pyplot as plt


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
        if i == 0:
            for j in range(1, x_i + 1):
                f[i + 1][j] = f[i][j] - c * 0.5 * (f[i][j + 1] - f[i][j - 1]) + d * (f[i][j + 1] - 2 * f[i][j] + f[i][j - 1])
                if j == x_i - 1:
                    f[i + 1][j + 2] = f[i + 1][j]
        else:
            for j in range(1, x_i + 1):
                f[i + 1][j] = f[i][j] + 1.5 * (- c * 0.5 * (f[i][j + 1] - f[i][j - 1]) + d * (f[i][j + 1] - 2 * f[i][j] + f[i][j - 1]))\
                              - 0.5 * (- c * 0.5 * (f[i - 1][j + 1] - f[i - 1][j - 1]) + d * (f[i - 1][j + 1] - 2 * f[i - 1][j] + f[i - 1][j - 1]))

                # f[i + 1][j] = f[i][j] - c * 0.5 * (f[i][j + 1] - f[i][j - 1]) + d * (f[i][j + 1] - 2 * f[i][j] + f[i][j - 1])
                if j == x_i - 1:
                    f[i + 1][j + 2] = f[i + 1][j]

    # graph
    fig, ax = plt.subplots()
    # p1 = plt.plot(x, f[0][:-1], marker="o", clip_on=False, label="t = 0")
    # p2 = plt.plot(x, f[100][:-1], marker="+", ms=8, clip_on=False, label="t = 100")
    # p3 = plt.plot(x, f[200][:-1], marker="s", clip_on=False, label="t = 200")
    # p4 = plt.plot(x, f[300][:-1], marker="D", clip_on=False, label="t = 300")
    # p5 = plt.plot(x, f[400][:-1], marker="x", clip_on=False, label="t = 400")
    # p6 = plt.plot(x, f[500][:-1], marker="^", clip_on=False, label="t = 500")
    p1 = plt.plot(x, f[0][:-1], label="t = 0")
    p2 = plt.plot(x, f[20][:-1], label=f"t = {int(dt * 20)}")
    p3 = plt.plot(x, f[40][:-1], label=f"t = {int(dt * 40)}")
    p4 = plt.plot(x, f[60][:-1], label=f"t = {int(dt * 60)}")
    p5 = plt.plot(x, f[80][:-1], label=f"t = {int(dt * 80)}")
    p6 = plt.plot(x, f[100][:-1], label=f"t = {int(dt * 100)}")
    plt.title("t: 2-order(AB method)")
    plt.legend(loc='upper right')
    plt.xlim(0, 5)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('f')
    plt.show()
