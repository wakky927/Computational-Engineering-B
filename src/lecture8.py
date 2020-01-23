import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation


def animate_f_true(nframe):
    plt.cla()
    plt.plot(x, f_true[nframe])
    plt.xlim(0, L)
    plt.ylim(-3 * A / c / n, 3 * A / c / n)
    plt.title(f"hogehoge (t = {nframe})")
    a, b = xx(0)
    plt.xticks(a, b)
    ax.set_xticks(np.linspace(0, L, 2 * n + 1), minor=True)
    plt.yticks([-2 * A / c / n, 0.0, 2 * A / c / n], [f"{round(-2 * A / c / n, 3)}", "0.0", f"{round(2 * A / c / n, 3)}"])
    ax.set_yticks(np.linspace(-3 * A / c / n, 3 * A / c / n, 3), minor=True)
    plt.xlabel('x')
    plt.ylabel('f')
    ax.grid(which="both")


def animate_g_true(nframe):
    plt.cla()
    plt.plot(x, g_true[nframe], color="orange")
    plt.xlim(0, L)
    plt.ylim(-5 * A / c / (2 * n - 1), 5 * A / c / (2 * n - 1))
    plt.title(f"hogehoge (t = {nframe})")
    a, b = xx(1)
    plt.xticks(a, b)
    ax.set_xticks(np.linspace(0, L, 2 * n), minor=True)
    plt.yticks([-4 * A / (2 * n - 1) / c, 0.0, 4 * A / (2 * n - 1) / c], [f"{round(-4 * A / (2 * n - 1) / c, 3)}", "0.0", f"{round(4 * A / (2 * n - 1) / c, 3)}"])
    ax.set_yticks(np.linspace(-5 * A / c / (2 * n - 1), 5 * A / c / (2 * n - 1), 3), minor=True)
    plt.xlabel('x')
    plt.ylabel('f')
    ax.grid(which="both")


def animate_f(nframe):
    plt.cla()
    plt.plot(x, f[nframe])
    plt.xlim(0, L)
    plt.ylim(-3 * A / c / n, 3 * A / c / n)
    plt.title(f"hogehoge (t = {nframe})")
    a, b = xx(0)
    plt.xticks(a, b)
    ax.set_xticks(np.linspace(0, L, 2 * n + 1), minor=True)
    plt.yticks([-2 * A / c / n, 0.0, 2 * A / c / n], [f"{round(-2 * A / c / n, 3)}", "0.0", f"{round(2 * A / c / n, 3)}"])
    ax.set_yticks(np.linspace(-3 * A / c / n, 3 * A / c / n, 3), minor=True)
    plt.xlabel('x')
    plt.ylabel('f')
    ax.grid(which="both")


def animate_g(nframe):
    plt.cla()
    plt.plot(x, g[nframe][:-1], color="orange")
    plt.xlim(0, L)
    plt.ylim(-3 * A / c / n, 3 * A / c / n)
    plt.title(f"hogehoge (t = {nframe})")
    a, b = xx(1)
    plt.xticks(a, b)
    ax.set_xticks(np.linspace(0, L, 2 * n), minor=True)
    plt.yticks([-4 * A / (2 * n - 1) / c, 0.0, 4 * A / (2 * n - 1) / c], [f"{round(-4 * A / (2 * n - 1) / c, 3)}", "0.0", f"{round(4 * A / (2 * n - 1) / c, 3)}"])
    ax.set_yticks(np.linspace(-5 * A / c / (2 * n - 1), 5 * A / c / (2 * n - 1), 3), minor=True)
    plt.xlabel('x')
    plt.ylabel('f')
    ax.grid(which="both")


def xx(s):
    if s == 0:
        x1 = [0]
        x2 = ["0.0"]
        for k in range(1, 2 * n + 1):
            x1.append(k * PI / n)
            if k == 1:
                x2.append(f"pi/{n}")

            elif k == n:
                x2.append("pi")

            elif k == 2 * n:
                x2.append("2pi")

            else:
                x2.append(f"{k}pi/{n}")

    else:
        x1 = np.zeros(2 * n)
        x2 = ["0.0"]
        for k in range(1, 2 * n):
            x1[k] = 2 * k * PI / (2 * n - 1)
            if k == 2 * n - 1:
                x2.append("2pi")

            else:
                x2.append(f"{int(2 * k)}pi/{int(2 * n - 1)}")

    return x1, x2


# setup
PI = np.pi
n = 3
L = 2 * PI
x_i = 100
t_i = 100
dt1 = L / t_i
dt2 = L / t_i
dx = L / x_i
c = 1
A = 1

x = np.zeros(x_i + 1)
for i in range(x_i + 1):
    x[i] = dx * i

f_true = np.zeros((t_i + 1, x_i + 1))  # exercise 1
g_true = np.zeros((t_i + 1, x_i + 1))  # exercise 2
for i in range(t_i + 1):
    for j in range(x_i + 1):
        f_true[i][j] = 2 * A / c / n * np.sin(n / 2 * dx * j) * np.sin(c * n / 2 * dt1 * i)
        g_true[i][j] = 4 * A / c / (2 * n - 1) * np.sin((2 * n - 1) / 4 * dx * j) * np.sin((2 * n - 1) * c / 4 * dt2 * i)

# f = np.zeros((t_i + 1, x_i + 1))

if __name__ == '__main__':
    print("start!")

    f = np.zeros((t_i + 1, x_i + 1))  # exercise 1, include I.C. f(x, 0) = 0
    for i in range(1, t_i + 1):
        if i == 1:  # I.C.
            for j in range(1, x_i + 1):  # B.C. f(0, t) = 0
                if j == x_i:
                    f[i][j] = 0  # B.C. f(2π, t) = 0

                else:
                    f[i][j] = A * np.sin(n * dx * j / 2) * dt1  # I.C. f'(x, 0) = A * sin(n * x / 2)

        else:
            for j in range(1, x_i + 1):  # B.C. f(0, t) = 0
                if j == x_i:
                    f[i][j] = 0  # B.C. f(2π, t) = 0

                else:
                    f[i][j] = 2 * f[i - 1][j] - f[i - 2][j] + c**2 * dt1**2 / dx**2 * (f[i - 1][j + 1] - 2 * f[i - 1][j] + f[i - 1][j - 1])

    g = np.zeros((t_i + 1, x_i + 2))  # exercise 2, include I.C. g(x, 0) = 0
    for i in range(1, t_i + 1):
        if i == 1:  # I.C.
            for j in range(1, x_i + 2):  # B.C. f(0, t) = 0
                if j == x_i - 1:  # B.C. f'(2π, t) = 0 | 2nd Order
                    g[i][j] = A * np.sin((2 * n - 1) * dx * j / 4) * dt2  # I.C. f'(x, 0) = A * sin((2 * n - 1) * x / 4)
                    g[i][j + 2] = g[i][j]

                else:
                    g[i][j] = A * np.sin((2 * n - 1) * dx * j / 4) * dt2  # I.C. f'(x, 0) = A * sin((2 * n - 1) * x / 4)

        else:
            for j in range(1, x_i + 1):  # B.C. f(0, t) = 0
                if j == x_i - 1:  # B.C. f'(2π, t) = 0 | 2nd Order
                    g[i][j] = 2 * g[i - 1][j] - g[i - 2][j] + c**2 * dt1**2 / dx**2 * (g[i - 1][j + 1] - 2 * g[i - 1][j] + g[i - 1][j - 1])
                    g[i][j + 2] = g[i][j]

                else:
                    g[i][j] = 2 * g[i - 1][j] - g[i - 2][j] + c ** 2 * dt1 ** 2 / dx ** 2 * (g[i - 1][j + 1] - 2 * g[i - 1][j] + g[i - 1][j - 1])

    # graph
    # fig, ax = plt.subplots()
    # p0 = plt.plot(x, f_true[0], label="t = 0")
    # p1 = plt.plot(x, f_true[int(t_i * 0.2)], label=f"t = {dt * t_i * 0.2}")
    # p2 = plt.plot(x, f_true[int(t_i * 0.4)], label=f"t = {dt * t_i * 0.4}")
    # p3 = plt.plot(x, f_true[int(t_i * 0.6)], label=f"t = {dt * t_i * 0.6}")
    # p4 = plt.plot(x, f_true[int(t_i * 0.8)], label=f"t = {dt * t_i * 0.8}")
    # p5 = plt.plot(x, f_true[int(t_i * 1.0)], label=f"t = {dt * t_i * 1.0}")
    # plt.title("hogehoge")
    # plt.legend(loc='upper right')
    # plt.xlim(0, 2 * PI)
    # plt.ylim(-0.5, 0.5)
    # plt.xlabel('x')
    # plt.ylabel('f')
    # plt.show()

    os.makedirs('../data/lecture8/gif', exist_ok=True)
    fig, ax = plt.subplots()
    # true
    anim_f_true = animation.FuncAnimation(fig, animate_f_true, frames=t_i)
    anim_f_true.save('../data/lecture8/gif/1_true.gif', writer='imagemagick', fps=24)
    anim_g_true = animation.FuncAnimation(fig, animate_g_true, frames=t_i)
    anim_g_true.save('../data/lecture8/gif/2_true.gif', writer='imagemagick', fps=24)
    # approximate
    anim_f = animation.FuncAnimation(fig, animate_f, frames=t_i)
    anim_f.save('../data/lecture8/gif/1_approximate.gif', writer='imagemagick', fps=24)
    anim_g = animation.FuncAnimation(fig, animate_g, frames=t_i)
    anim_g.save('../data/lecture8/gif/2_approximate.gif', writer='imagemagick', fps=24)

    print("fin")
