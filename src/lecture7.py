import copy

import numpy as np
from matplotlib import pyplot as plt


u = 5.0
D = 1.0
n = 10
L = 2.0
dx = L / n
Pe = u * dx / D
omega = 1.6
EPS = 0.01
i_max = 1000

aw = u / dx + D / dx**2
ap = - u / dx - 2 * D / dx**2
ae = D / dx**2
# aw = 3
# ap = -4
# ae = 1

a = np.full(n + 1, aw)
b = np.full(n + 1, ap)
c = np.full(n + 1, ae)
a[-1] = c[0] = 0
b[0] = b[-1] = 1

f = np.zeros(n + 1)
x = np.zeros(n + 1)
for i in range(n + 1):
    x[i] = dx * i
    f[i] = (np.exp(Pe * dx * i / L) - 1) / (np.exp(Pe) - 1)

f_TDMA = np.zeros(n + 1)
f_SOR = np.zeros(n + 1)

b_TDMA = np.zeros(n + 1)
b_SOR = np.zeros(n + 1)
b_TDMA[-1] = 1
b_SOR[-1] = 1

if __name__ == '__main__':
    # TDMA method
    gam = np.zeros(n)
    bet = b[0]
    f_TDMA[0] = b_TDMA[0] / bet
    for i in range(1, n + 1):
        gam[i - 1] = c[i - 1] / bet
        bet = b[i] - a[i] * gam[i - 1]
        f_TDMA[i] = (b_TDMA[i] - a[i] * f_TDMA[i - 1]) / bet

    for j in reversed(range(1, n)):
        f_TDMA[j] -= gam[j] * f_TDMA[j + 1]

    # SOR method
    for i in range(i_max):
        M = 0
        f_tmp = copy.copy(f_SOR)
        for j in range(n + 1):
            if j == 0 or j == n:
                f_SOR[j] = f_tmp[j] * (1 - omega) + b_SOR[j] * omega
                r1 = np.abs((f_SOR[j] - f_tmp[j]) / f_SOR[j])
                if r1 < EPS:
                    M += 1

            else:
                f_SOR[j] = f_tmp[j] * (1 - omega) + (b_SOR[j] - c[j] * f_tmp[j + 1] - a[j] * f_SOR[j - 1]) / b[j] * omega
                r2 = np.abs((f_SOR[j] - f_tmp[j]) / f_SOR[j])
                if r2 < EPS:
                    M += 1

        if M == n:
            break

    print("fin")

    # graph1
    fig1, ax1 = plt.subplots()
    p1 = plt.plot(x, f, label="true")
    p2 = plt.plot(x, f_TDMA, label="TDMA")
    plt.title("TDMA method")
    plt.legend(loc='upper left')
    plt.xlim(0, 2)
    plt.ylim(0, 1)
    plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0], ["0.0", "0.lecture10", "0.4", "0.6", "0.8", "1.0", "1.lecture10", "1.4", "1.6", "1.8", "lecture10.0"])
    ax1.set_xticks(np.linspace(0, 2, 11), minor=True)
    plt.yticks([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"])
    ax1.set_yticks(np.linspace(0, 1, 3), minor=True)
    plt.xlabel('x')
    plt.ylabel('f')
    plt.show()

    # graph2
    fig2, ax2 = plt.subplots()
    q1 = plt.plot(x, f, label="true")
    q2 = plt.plot(x, f_SOR, label="SOR", color="green")
    plt.title("TDMA method")
    plt.legend(loc='upper left')
    plt.xlim(0, 2)
    plt.ylim(0, 1)
    plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0], ["0.0", "0.lecture10", "0.4", "0.6", "0.8", "1.0", "1.lecture10", "1.4", "1.6", "1.8", "lecture10.0"])
    ax2.set_xticks(np.linspace(0, 2, 11), minor=True)
    plt.yticks([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"])
    ax2.set_yticks(np.linspace(0, 1, 3), minor=True)
    plt.xlabel('x')
    plt.ylabel('f')
    plt.show()
