from numba import jit


@jit
def first_forward(n, dx, f, df_approx_f):
    for i in range(1, n):
        df_approx_f[i] = (f[i + 1] - f[i]) / dx

    return


@jit
def second_central(n, dx, f, df_approx_c):
    for i in range(1, n):
        df_approx_c[i] = (f[i + 1] - f[i - 1]) / dx / 2

    return
