from numba import jit


''' Lecture1 '''
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


''' Lecture2 '''
@jit
def lec2_first_func(dx, f):
    f1_approx = (0.5 * f[0] + 0.5 * f[1]) * dx

    return f1_approx


@jit
def lec2_second_func(dx, f):
    f2_approx = (1 * f[0]) * dx

    return f2_approx


''' Lecture3 '''
@jit
def lec3_euler(n, dt, f_approx):
    f1_approx = f_approx[1]

    for i in range(2, n):
        f_approx[i] = f1_approx + dt * (1.0 * f1_approx)
        f1_approx = f_approx[i]

    return


@jit
def lec3_ab(n, dt, f_approx):
    f1_approx = f_approx[1]
    f2_approx = f_approx[0]

    for i in range(2, n):
        f_approx[i] = f1_approx + dt * (1.5 * f1_approx - 0.5 * f2_approx)
        f2_approx = f1_approx
        f1_approx = f_approx[i]

    return


@jit
def lec3_tra(n, dt, f_approx):
    f1_approx = f_approx[1]

    for i in range(2, n):
        f_approx[i] = f1_approx + dt * (0.5 * f1_approx + 0.5 * f_approx[i])
        f1_approx = f_approx[i]

    return
