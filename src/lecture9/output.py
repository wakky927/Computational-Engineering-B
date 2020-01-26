import os

import numpy as np


def grid(xp, yp):
    os.makedirs('../../data/lecture9', exist_ok=True)
    grid_mat = np.stack([xp, yp])
    np.savetxt(
        '../../data/lecture9/grid.csv', grid_mat, delimiter=',', fmt='%.10f')

    return


def solution(p, u, v, m, n):
    print(f"\nvelocity u")
    print(u)

    print(f"\nvelocity v")
    print(v)

    print(f"\npressure")
    print(p)

    return


def divergent(p, u, v, dx, dy, m, n):
    md = 101
    nd = 101

    div = np.zeros((md, nd))

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # div[i][j] = (u[i + 1][j] - u[i - 1][j]) / dx / 2\
            #             + (v[i][j + 1] - v[i][j - 1]) / dy / 2

            div[i][j] = (u[i][j] - u[i - 1][j]) / dx\
                        + (v[i][j] - v[i][j - 1]) / dy

    os.makedirs('../../data/lecture9', exist_ok=True)
    np.savetxt(
        '../../data/lecture9/divergent.csv', div, delimiter=',', fmt='%.10f')

    return


def solution_post(p, u, v, m, n):
    md = 101
    nd = 101

    u_cnt = np.zeros((md, nd))
    v_cnt = np.zeros((md, nd))

    # interpolation at p-center grid
    # for i in range(1, m + 1):
    #     for j in range(1, n + 1):
    #         u_cnt[i][j] = u[i][j]
    #         v_cnt[i][j] = v[i][j]
    #
    # for j in range(1, n + 1):
    #     u_cnt[0][j] = u_cnt[1][j]
    #     v_cnt[0][j] = v_cnt[1][j]
    #     u_cnt[m + 1][j] = u_cnt[m][j]
    #     v_cnt[m + 1][j] = v_cnt[m][j]
    #
    # for i in range(m + 2):
    #     u_cnt[i][0] = u_cnt[i][1]
    #     v_cnt[i][0] = v_cnt[i][1]
    #     u_cnt[i][n + 1] = u_cnt[i][n]
    #     v_cnt[i][n + 1] = v_cnt[i][n]

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            u_cnt[i][j] = 0.5 * (u[i][j] + u[i - 1][j])
            v_cnt[i][j] = 0.5 * (v[i][j] + v[i][j - 1])

    for j in range(1, n + 1):
        u_cnt[0][j] = u[0][j]
        v_cnt[0][j] = 0.5 * (v[0][j] + v[0][j - 1])
        u_cnt[m + 1][j] = 0.5 * (u[m + 1][j] + u[m][j])
        v_cnt[m + 1][j] = 0.5 * (v[m + 1][j] + v[m + 1][j - 1])

    for i in range(m + 2):
        u_cnt[i][0] = 0.5 * (u[i][0] + u[i - 1][0])
        v_cnt[i][0] = v[i][0]
        u_cnt[i][n + 1] = 0.5 * (u[i][n + 1] + u[i - 1][n + 1])
        v_cnt[i][n + 1] = v[i][n]

    os.makedirs('../../data/lecture9', exist_ok=True)
    np.savetxt(
        '../../data/lecture9/velocity_u.csv', u_cnt, delimiter=',',
        fmt='%.10f')
    np.savetxt(
        '../../data/lecture9/velocity_v.csv', v_cnt, delimiter=',',
        fmt='%.10f')
    np.savetxt(
        '../../data/lecture9/velocity_u.csv', p, delimiter=',', fmt='%.10f')


    return


def paraview(p, xp, yp, m, n, u, v):
    os.makedirs('../../data/lecture9', exist_ok=True)
    path_w = '../../data/lecture9/output_paraview.vtk'

    with open(path_w, mode='w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("2D flow\n")
        f.write("ASCII\n")

        f.write("DATASET STRUCTURED_GRID\n")
        f.write(f"DIMENSIONS {m} {n} 1\n")

        f.write(f"POINTS {m * n} float\n")
        for j in range(1, n + 1):
            for i in range(1, m + 1):
                f.write(f"{round(xp[i], 4)} {round(yp[j], 4)} 0.0000\n")

        f.write(f"POINT_DATA {m * n}\n")

        # velocity vector
        f.write(f"VECTORS velocity float\n")
        for j in range(1, n + 1):
            for i in range(1, m + 1):
                f.write(f"{round(u[i][j], 4)} {round(v[i][j], 4)} 0.0000\n")

        # pressure
        f.write(f"SCALARS pressure float\n")
        f.write(f"LOOKUP_TABLE default\n")
        for j in range(1, n + 1):
            for i in range(1, m + 1):
                f.write(f"{round(p[i][j], 4)}\n")

    return
