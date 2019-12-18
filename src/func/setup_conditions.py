import numpy as np
from numba import jit


@jit
def lec1(n, dx, x, f, df_exact):
    for i in range(n + 1):
        x[i] = dx * i
        f[i] = np.exp(x[i])
        df_exact[i] = f[i]

    return


@jit
def lec2(n, dx, x, f):
    for i in range(n + 1):
        x[i] = dx * i
        f[i] = np.exp(x[i])

    f1_exact = np.exp(dx) - np.exp(0)
    f2_exact = np.exp(dx / 2) - np.exp(-dx / 2)

    return f1_exact, f2_exact


@jit
def lec3(n, dt, t, f_exact):
    for i in range(n + 1):
        t[i] = dt * i
        f_exact[i] = np.exp(t[i])

    return
