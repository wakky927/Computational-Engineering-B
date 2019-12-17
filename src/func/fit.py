import numpy as np
from numba import jit


@jit
def xp(x):
    return np.linspace(np.min(x), np.max(x), 100)


@jit
def loglog(x, a, b):
    return a * x**b
