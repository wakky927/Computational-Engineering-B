import numpy as np
from numba import jit


@jit
def SOR(md, p, ap, ae, aw, bb, m, p_exact, relax_factor):
    iter_n, error1, error2, error3 = 0, 0, 0, 0  # initialize parameters
    p_old = np.zeros(md + 1)

    ''' SOR algorithm '''
    iter_max = 50  # SOR max iteration steps

    for iter_n in range(1, iter_max + 1):
        error1 = 0  # error reset
        error2 = 0
        error3 = 0

        for i in range(1, m + 1):  # step increase
            p_old[i] = p[i]

        for i in range(1, m + 1):
            p_diff = -p_old[i] + (bb[i] - ae[i] * p_old[i + 1] - aw[i] * p[i - 1]) / ap[i]
            p[i] = p_old[i] + p_diff * relax_factor

            error1 = max(error1, abs(p[i] - p_old[i]))
            error2 = max(error2, abs(p[i] - p_exact[i]))
            error3 = max(error3, abs(p_diff))

    return iter_n, error1, error2, error3
