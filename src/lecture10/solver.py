import numpy as np
from numba import jit

import condition


@jit
def solve_matrix(p, Ap, Ae, Aw, An, As, bb, m, n):
    md = 202
    nd = 202

    p_old = np.zeros((md, nd))

    ''' SOR algorithm '''
    iter_max = 300  # SOR max iteration steps
    relax_factor = 1.8  # SOR relaxation factor

    for iter_i in range(1, iter_max + 1):
        error = 0.0

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                p_old[i][j] = p[i][j]

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                p[i][j] = (bb[i][j] - Ae[i][j] * p_old[i + 1][j] - Aw[i][j]
                           * p[i - 1][j] - An[i][j] * p_old[i][j + 1]
                           - As[i][j] * p[i][j - 1]) / Ap[i][j]\
                          * relax_factor + p_old[i][j] * (1 - relax_factor)
                a = np.abs(p[i][j] - p_old[i][j])
                e = max(error, a)
                error = e

    # print(f"iteration no. {iter_i} -- error = {error}")  # jit cannot print

    return


@jit
def solve_p(p, u, v, u_old, v_old, nue, density, dx, dy, dt, m, n, xp, yp,
            height, length):
    md = 202
    nd = 202

    Ap = np.zeros((md, nd))
    Ae = np.zeros((md, nd))
    Aw = np.zeros((md, nd))
    An = np.zeros((md, nd))
    As = np.zeros((md, nd))
    bb = np.zeros((md, nd))

    # u_stg = 0.0
    # v_stg = 0.0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            ''' velocity u '''
            # convection_x (1st upwind scheme)
            u[i][j] = u_old[i][j]\
                      - dt * max(u_old[i][j], 0.0)\
                      * (u_old[i][j] - u_old[i - 1][j]) / dx\
                      - dt * min(u_old[i][j], 0.0)\
                      * (u_old[i + 1][j] - u_old[i][j]) / dx

            # convection_y
            # v_stg = 0.25 * (v_old[i][j] + v_old[i + 1][j] + v_old[i][j - 1]
            #                 + v_old[i + 1][j - 1])  # Staggered grid
            u[i][j] = u[i][j]\
                      - dt * max(v_old[i][j], 0.0)\
                      * (u_old[i][j] - u_old[i][j - 1]) / dy\
                      - dt * min(v_old[i][j], 0.0)\
                      * (u_old[i][j + 1] - u_old[i][j]) / dy

            # diffusion_x
            u[i][j] = u[i][j]\
                      + dt * nue * (u_old[i + 1][j] - 2 * u_old[i][j]
                                    + u_old[i - 1][j]) / dx**2

            # diffusion_y
            u[i][j] = u[i][j] \
                      + dt * nue * (u_old[i][j + 1] - 2 * u_old[i][j]
                                    + u_old[i][j - 1]) / dy**2

            ''' velocity v '''
            # convection_x (1st upwind scheme)
            # u_stg = 0.25 * (u_old[i][j] + u_old[i - 1][j] + u_old[i][j + 1]
            #                 + u_old[i - 1][j + 1])  # Staggered grid
            v[i][j] = v_old[i][j] \
                      - dt * max(u_old[i][j], 0.0) \
                      * (v_old[i][j] - v_old[i - 1][j]) / dx \
                      - dt * min(u_old[i][j], 0.0) \
                      * (v_old[i + 1][j] - v_old[i][j]) / dx

            # convection_y
            v[i][j] = v[i][j] \
                      - dt * max(v_old[i][j], 0.0) \
                      * (v_old[i][j] - v_old[i][j - 1]) / dy \
                      - dt * min(v_old[i][j], 0.0) \
                      * (v_old[i][j + 1] - v_old[i][j]) / dy

            # diffusion_x
            v[i][j] = v[i][j] \
                      + dt * nue * (v_old[i + 1][j] - 2 * v_old[i][j]
                                    + v_old[i - 1][j]) / dx**2

            # diffusion_y
            u[i][j] = u[i][j] \
                      + dt * nue * (v_old[i][j + 1] - 2 * v_old[i][j]
                                    + v_old[i][j - 1]) / dy**2

    ''' matrix solution '''
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            Ae[i][j] = dt / density / dx**2
            Aw[i][j] = dt / density / dx**2
            An[i][j] = dt / density / dy**2
            As[i][j] = dt / density / dy**2
            Ap[i][j] = - Ae[i][j] - Aw[i][j] - An[i][j] - As[i][j]

            # bb[i][j] = (u[i][j] - u[i - 1][j]) / dx\
            #            + (v[i][j] - v[i][j - 1]) / dy
            bb[i][j] = (u[i + 1][j] - u[i - 1][j]) / dx / 2 \
                       + (v[i][j + 1] - v[i][j - 1]) / dy / 2

    condition.matrix_c(p, Ap, Ae, Aw, An, As, bb, m, n, xp, yp, height, length)
    solve_matrix(p, Ap, Ae, Aw, An, As, bb, m, n)

    return


@jit
def solve_u(p, u, density, dx, dt, m, n):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # convection_x (1st upwind scheme) -> already calculated in solve_p
            # convection_y -> already calculated in solve_p
            # diffusion_x -> already calculated in solve_p
            # diffusion_y -> already calculated in solve_p
            # pressure
            # u[i][j] = u[i][j]\
            #           - dt / density * (p[i + 1][j] - p[i][j]) / dx
            u[i][j] = u[i][j] \
                      - dt / density * (p[i + 1][j] - p[i - 1][j]) / dx / 2


@jit
def solve_v(p, v, density, dy, dt, m, n):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # convection_x (1st upwind scheme) -> already calculated in solve_p
            # convection_y -> already calculated in solve_p
            # diffusion_x -> already calculated in solve_p
            # diffusion_y -> already calculated in solve_p
            # pressure
            # v[i][j] = v[i][j] \
            #           - dt / density * (p[i][j + 1] - p[i][j]) / dy
            v[i][j] = v[i][j] \
                      - dt / density * (p[i][j + 1] - p[i][j - 1]) / dy / 2
