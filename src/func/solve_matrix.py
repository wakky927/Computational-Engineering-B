import numpy as np
from numba import jit


@jit
def SOR1(md, p, ap, ae, aw, bb, m, p_exact, relax_factor):
    eps = 1e-15  # convergence criterion
    error1, error2, error3 = 0, 0, 0  # initialize parameters
    p_old = np.zeros(md + 1)

    ''' SOR algorithm '''
    iter_n = 1

    while True:
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

        if error1 < eps:
            break

        elif iter_n > 500:
            break

        else:
            iter_n += 1

    return iter_n, error1, error2, error3


@jit
def SOR3(md, p, ap, ae, aw, bb, m, p_exact, relax_factor):
    eps = 1e-15  # convergence criterion
    error1, error2, error3 = 0, 0, 0  # initialize parameters
    p_old = np.zeros(md + 1)

    ''' SOR algorithm '''
    iter_n = 1

    while True:
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

        if error3 < eps:
            break

        elif iter_n > 500:
            break

        else:
            iter_n += 1

    return iter_n, error1, error2, error3
