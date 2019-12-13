import numpy as np
from numba import jit


@jit
def lec1(n, dx, x, f, df_exact):
    for i in range(n + 1):
        x[i] = dx * i
        f[i] = np.exp(x[i])
        df_exact[i] = f[i]

    return
