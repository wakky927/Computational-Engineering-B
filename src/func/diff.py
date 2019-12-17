from numba import jit


@jit
def lec1_first_forward(n, dx, f, df_approx_f):
    for i in range(1, n):
        df_approx_f[i] = (f[i + 1] - f[i]) / dx

    return


@jit
def lec1_second_central(n, dx, f, df_approx_c):
    for i in range(1, n):
        df_approx_c[i] = (f[i + 1] - f[i - 1]) / dx / 2

    return


@jit
def lec2_first_func(dx, f):
    f1_approx = (0.5 * f[0] + 0.5 * f[1]) * dx

    return f1_approx


@jit
def lec2_second_func(dx, f):
    f2_approx = (1 * f[0]) * dx

    return f2_approx
