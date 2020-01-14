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
PE1 = 1
PE2 = 2.0
PE3 = 3.0
D = 1.0e-9
L = 10.0
dx = 0.05
u1 = PE1 * D / dx
u2 = PE2 * D / dx
u3 = PE3 * D / dx
dt = dt2 = 625000
dt1 = dt2 * PE2 / PE1
dt3 = dt2 * PE2 / PE3
# c1 = u1 * dt1 / dx
# c2 = u2 * dt2 / dx
# c3 = u3 * dt3 / dx
c = c1 = c2 = c3 = 0.5
d1 = D * dt1 / dx / dx  # 0.5
d2 = D * dt2 / dx / dx  # 0.25
d3 = D * dt3 / dx / dx  # 0.167

t_i = 1501
x_i = int(L / dx)
x = np.zeros(x_i + 1)
for i in range(x_i + 1):
    x[i] = dx * i
f = np.zeros((t_i + 1, x_i + 2))
f1 = np.zeros((t_i + 1, x_i + 2))
f2 = np.zeros((t_i + 1, x_i + 2))
f3 = np.zeros((t_i + 1, x_i + 2))

# I.C.
for i in range(int(1.0 / dx) + 1):
    f[0][i] = 0.5 * (1 - np.cos(2 * np.pi * x[i]))
    f1[0][i] = 0.5 * (1 - np.cos(2 * np.pi * x[i]))
    f2[0][i] = 0.5 * (1 - np.cos(2 * np.pi * x[i]))
    f3[0][i] = 0.5 * (1 - np.cos(2 * np.pi * x[i]))

if __name__ == '__main__':
    for i in range(t_i):
        for j in range(1, x_i + 1):
            f[i + 1][j] = f[i][j] - c * (f[i][j] - f[i][j - 1])  # Convention eq.
            f1[i + 1][j] = f1[i][j] - c1 * 0.5 * (f1[i][j + 1] - f1[i][j - 1]) + d1 * (f1[i][j + 1] - 2 * f1[i][j] + f1[i][j - 1])
            f2[i + 1][j] = f2[i][j] - c2 * 0.5 * (f2[i][j + 1] - f2[i][j - 1]) + d2 * (f2[i][j + 1] - 2 * f2[i][j] + f2[i][j - 1])
            f3[i + 1][j] = f3[i][j] - c3 * 0.5 * (f3[i][j + 1] - f3[i][j - 1]) + d3 * (f3[i][j + 1] - 2 * f3[i][j] + f3[i][j - 1])
            if j == x_i - 1:  # Norman condition
                f[i + 1][j + 2] = f[i + 1][j]
                f1[i + 1][j + 2] = f1[i + 1][j]
                f2[i + 1][j + 2] = f2[i + 1][j]
                f3[i + 1][j + 2] = f3[i + 1][j]

    # graph1
    fig1, ax1 = plt.subplots()
    p1_1 = plt.plot(x, f1[0][:-1], label="t = 0")
    p1_2 = plt.plot(x, f1[20][:-1], label=f"t = {round(dt1 * 20 / 10000000, 2)}E+7")
    p1_3 = plt.plot(x, f1[40][:-1], label=f"t = {round(dt1 * 40 / 10000000, 2)}E+7")
    p1_4 = plt.plot(x, f1[60][:-1], label=f"t = {round(dt1 * 60 / 10000000, 2)}E+7")
    p1_5 = plt.plot(x, f1[80][:-1], label=f"t = {round(dt1 * 80 / 100000000, 2)}E+8")
    p1_6 = plt.plot(x, f1[100][:-1], label=f"t = {round(dt1 * 100 / 100000000, 2)}E+8")
    q1_1 = plt.plot(x, f[0][:-1], linestyle="dashed", color="blue", label="t = 0")
    q1_2 = plt.plot(x, f[40][:-1], linestyle="dashed", color="orange", label=f"t = {round(dt * 40 / 10000000, 2)}E+7")
    q1_3 = plt.plot(x, f[80][:-1], linestyle="dashed", color="green", label=f"t = {round(dt * 80 / 10000000, 2)}E+7")
    q1_4 = plt.plot(x, f[120][:-1], linestyle="dashed", color="red", label=f"t = {round(dt * 120 / 10000000, 2)}E+7")
    q1_5 = plt.plot(x, f[160][:-1], linestyle="dashed", color="purple", label=f"t = {round(dt * 160 / 100000000, 2)}E+8")
    q1_6 = plt.plot(x, f[200][:-1], linestyle="dashed", color="brown", label=f"t = {round(dt * 200 / 100000000, 2)}E+8")
    plt.title(f"Convective diffusion eq. (Pe = {PE1}) & Convention eq.", fontsize=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(right=0.7)
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('f')
    plt.show()

    # graph2
    fig2, ax2 = plt.subplots()
    p2_1 = plt.plot(x, f2[0][:-1], label="t = 0")
    p2_2 = plt.plot(x, f2[40][:-1], label=f"t = {round(dt2 * 40 / 10000000, 2)}E+7")
    p2_3 = plt.plot(x, f2[80][:-1], label=f"t = {round(dt2 * 80 / 10000000, 2)}E+7")
    p2_4 = plt.plot(x, f2[120][:-1], label=f"t = {round(dt2 * 120 / 10000000, 2)}E+7")
    p2_5 = plt.plot(x, f2[160][:-1], label=f"t = {round(dt2 * 160 / 100000000, 2)}E+8")
    p2_6 = plt.plot(x, f2[200][:-1], label=f"t = {round(dt2 * 200 / 100000000, 2)}E+8")
    q2_1 = plt.plot(x, f[0][:-1], linestyle="dashed", color="blue", label="t = 0")
    q2_2 = plt.plot(x, f[40][:-1], linestyle="dashed", color="orange", label=f"t = {round(dt * 40 / 10000000, 2)}E+7")
    q2_3 = plt.plot(x, f[80][:-1], linestyle="dashed", color="green", label=f"t = {round(dt * 80 / 10000000, 2)}E+7")
    q2_4 = plt.plot(x, f[120][:-1], linestyle="dashed", color="red", label=f"t = {round(dt * 120 / 10000000, 2)}E+7")
    q2_5 = plt.plot(x, f[160][:-1], linestyle="dashed", color="purple", label=f"t = {round(dt * 160 / 100000000, 2)}E+8")
    q2_6 = plt.plot(x, f[200][:-1], linestyle="dashed", color="brown", label=f"t = {round(dt * 200 / 100000000, 2)}E+8")
    plt.title(f"Convective diffusion eq. (Pe = {PE2}) & Convention eq.", fontsize=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(right=0.7)
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('f')
    plt.show()

    # graph3
    fig3, ax3 = plt.subplots()
    p3_1 = plt.plot(x, f3[0][:-1], label="t = 0")
    p3_2 = plt.plot(x, f3[60][:-1], label=f"t = {round(dt3 * 60 / 10000000, 2)}E+7")
    p3_3 = plt.plot(x, f3[120][:-1], label=f"t = {round(dt3 * 120 / 10000000, 2)}E+7")
    p3_4 = plt.plot(x, f3[180][:-1], label=f"t = {round(dt3 * 180 / 10000000, 2)}E+7")
    p3_5 = plt.plot(x, f3[240][:-1], label=f"t = {round(dt3 * 240 / 100000000, 2)}E+8")
    p3_6 = plt.plot(x, f3[300][:-1], label=f"t = {round(dt3 * 300 / 100000000, 2)}E+8")
    q3_1 = plt.plot(x, f[0][:-1], linestyle="dashed", color="blue", label="t = 0")
    q3_2 = plt.plot(x, f[40][:-1], linestyle="dashed", color="orange", label=f"t = {round(dt * 40 / 10000000, 2)}E+7")
    q3_3 = plt.plot(x, f[80][:-1], linestyle="dashed", color="green", label=f"t = {round(dt * 80 / 10000000, 2)}E+7")
    q3_4 = plt.plot(x, f[120][:-1], linestyle="dashed", color="red", label=f"t = {round(dt * 120 / 10000000, 2)}E+7")
    q3_5 = plt.plot(x, f[160][:-1], linestyle="dashed", color="purple", label=f"t = {round(dt * 160 / 100000000, 2)}E+8")
    q3_6 = plt.plot(x, f[200][:-1], linestyle="dashed", color="brown", label=f"t = {round(dt * 200 / 100000000, 2)}E+8")
    plt.title(f"Convective diffusion eq. (Pe = {PE3}) & Convention eq.", fontsize=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(right=0.7)
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('f')
    # plt.show()

    fig, ax = plt.subplots()
    anim = animation.FuncAnimation(fig, animate, frames=500)
    os.makedirs('../../data/lecture6/images', exist_ok=True)
    anim.save('convection_diffusion_6_3.gif', writer='imagemagick', fps=30)
