import os

import numpy as np
from matplotlib import pyplot as plt


def grid(xp, yp, m, n):
    os.makedirs(f'../../data/lecture9/{m}', exist_ok=True)
    grid_mat = np.stack([xp[1:m+1], yp[1:n+1]])
    np.savetxt(
        f'../../data/lecture9/{m}/grid.csv', grid_mat, delimiter=',',
        fmt='%.10f')

    return


def solution(p, u, v):
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
            div[i][j] = (u[i + 1][j] - u[i - 1][j]) / dx / 2\
                        + (v[i][j + 1] - v[i][j - 1]) / dy / 2

            # div[i][j] = (u[i][j] - u[i - 1][j]) / dx\
            #             + (v[i][j] - v[i][j - 1]) / dy

    os.makedirs(f'../../data/lecture9/{m}', exist_ok=True)
    np.savetxt(
        f'../../data/lecture9/{m}/divergent.csv', div, delimiter=',',
        fmt='%.10f')

    return


def solution_post(p, u, v, m, n, yp, height, inlet_velocity):
    md = 101
    nd = 101

    u_cnt = np.zeros((md, nd))
    v_cnt = np.zeros((md, nd))
    uu = np.zeros((m + 2, n + 2))
    vv = np.zeros((m + 2, n + 2))

    u_ans = np.zeros(m + 2)
    H = 2 * height
    for i in range(1, m + 2):
        u_ans[i] = - 6 / H**2 * inlet_velocity * (yp[i] - H) * yp[i]

    # interpolation at p-center grid
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            u_cnt[i][j] = u[i][j]
            v_cnt[i][j] = v[i][j]

    for j in range(1, n + 1):
        u_cnt[0][j] = u_cnt[1][j]
        v_cnt[0][j] = v_cnt[1][j]
        u_cnt[m + 1][j] = u_cnt[m][j]
        v_cnt[m + 1][j] = v_cnt[m][j]

    for i in range(m + 2):
        u_cnt[i][0] = u_cnt[i][1]
        v_cnt[i][0] = v_cnt[i][1]
        u_cnt[i][n + 1] = u_cnt[i][n]
        v_cnt[i][n + 1] = v_cnt[i][n]

    # for i in range(1, m + 1):
    #     for j in range(1, n + 1):
    #         u_cnt[i][j] = 0.5 * (u[i][j] + u[i - 1][j])
    #         v_cnt[i][j] = 0.5 * (v[i][j] + v[i][j - 1])
    #
    # for j in range(1, n + 1):
    #     u_cnt[0][j] = u[0][j]
    #     v_cnt[0][j] = 0.5 * (v[0][j] + v[0][j - 1])
    #     u_cnt[m + 1][j] = 0.5 * (u[m + 1][j] + u[m][j])
    #     v_cnt[m + 1][j] = 0.5 * (v[m + 1][j] + v[m + 1][j - 1])
    #
    # for i in range(m + lecture10):
    #     u_cnt[i][0] = 0.5 * (u[i][0] + u[i - 1][0])
    #     v_cnt[i][0] = v[i][0]
    #     u_cnt[i][n + 1] = 0.5 * (u[i][n + 1] + u[i - 1][n + 1])
    #     v_cnt[i][n + 1] = v[i][n]

    for i in range(1, m + 2):
        for j in range(1, n + 2):
            uu[i][j] = u_cnt[j][i]
            vv[i][j] = v_cnt[j][i]

    e = abs(u_ans[2:m + 1] - uu[2:m + 1, n]) / uu[2:m + 1, n] * 100

    os.makedirs(f'../../data/lecture9/{m}', exist_ok=True)
    np.savetxt(
        f'../../data/lecture9/{m}/velocity_u.csv', uu,
        delimiter=',',
        fmt='%.10f')
    np.savetxt(
        f'../../data/lecture9/{m}/velocity_v.csv', vv,
        delimiter=',',
        fmt='%.10f')
    np.savetxt(
        f'../../data/lecture9/{m}/pressure.csv', p,
        delimiter=',',
        fmt='%.10f')

    np.savetxt(f'../../data/lecture9/e{m}.csv', e, delimiter=',')

    fig1, ax1 = plt.subplots(figsize=(8, 6))
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 10
    ans = ax1.plot(yp[1:m + 1], u_ans[1:m + 1], label='true', color='blue')
    app = ax1.plot(yp[1:m + 1], uu[1:m + 1, n], label='approximate',
                   color='orange')
    ax2 = ax1.twinx()
    error = ax2.plot(yp[2:m + 1], e, label='error', linestyle="dashed",
                     color='green')
    plt.title(f"n = {n}, m = {m}")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc='lower right')
    ax1.set_xlabel('y [m]')
    ax1.set_ylabel('velocity u [m/s]')
    ax2.set_ylabel('error [%]')
    plt.show()
    fig1.savefig(f'../../data/lecture9/{m}/error.png')

    return


def paraview(p, xp, yp, m, n, u, v):
    os.makedirs(f'../../data/lecture9/{m}', exist_ok=True)
    path_w = f'../../data/lecture9/{m}/output_paraview.vtk'

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
