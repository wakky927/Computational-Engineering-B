import numpy as np
from numba import jit


''' Lecture1 '''
@jit
def lec1(n, dx, x, f, df_exact):
    for i in range(n + 1):
        x[i] = dx * i
        f[i] = np.exp(x[i])
        df_exact[i] = f[i]

    return


''' Lecture2 '''
@jit
def lec2(n, dx, x, f):
    for i in range(n + 1):
        x[i] = dx * i
        f[i] = np.exp(x[i])

    f1_exact = np.exp(dx) - np.exp(0)
    f2_exact = np.exp(dx / 2) - np.exp(-dx / 2)

    return f1_exact, f2_exact


''' Lecture3 '''
@jit
def lec3(n, dt, t, f_exact):
    for i in range(n + 1):
        t[i] = dt * i
        f_exact[i] = np.exp(t[i])

    return


''' Lecture4 '''
@jit
def lec4(m, dl, dx, c, x, ap, aw, ae, bb):
    for i in range(1, m + 1):
        x[i] = dl * (i - 1) / (m - 1)
        ap[i] = - 2
        aw[i] = 1.0
        ae[i] = 1.0
        bb[i] = c * dx**2

    return
