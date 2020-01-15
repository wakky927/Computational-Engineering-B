import copy

import numpy as np

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

a = np.full(n + 1, aw)
b = np.full(n + 1, ap)
c = np.full(n + 1, ae)
a[-1] = c[0] = 0
b[0] = b[-1] = 1

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
