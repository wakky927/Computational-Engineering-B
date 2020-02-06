from numba import jit


def physical_c():
    # nue = 1.0e-4  # [m^2/s] : water
    # density = 1000.0  # [kg/m^3] : water
    nue = 1.0e-3  # [m^2/s] : glycerin
    density = 1261.0  # [kg/m^3] : glycerin

    length = 0.1  # [m]
    height = 0.1  # [m]
    time = 100.0  # [s]
    inlet_velocity = 0.05  # [m/s]
    outlet_pressure = 0.0  # [N/m^lecture10]

    Reynolds_no = height * inlet_velocity / nue

    print(f"\nnue = {nue}")
    print(f"density = {density}")
    print(f"length = {length}")
    print(f"height = {height}")
    print(f"time = {time}")
    print(f"inlet_velocity = {inlet_velocity}")
    print(f"outlet_pressure = {outlet_pressure}")
    print(f"Reynolds_no = {Reynolds_no}")

    return nue, density, length, height, time, inlet_velocity, outlet_pressure


def grid_c(xp, yp, nue, length, height, time, inlet_velocity):
    mm = 10  # on body to x
    nn = 10  # on body to y
    m1 = 20  # front body to x
    m2 = 70  # behind body to x
    n1 = 30  # below body to y
    n2 = 30  # above body to y
    # mm = 5  # on body to x
    # nn = 5  # on body to y
    # m1 = 20  # front body to x
    # m2 = 175  # behind body to x
    # n1 = 1  # below body to y
    # n2 = 94  # above body to y
    m = mm + m1 + m2  # whole size to x (202)
    n = nn + n1 + n2  # whole size to y (102)
    istep_max = 5000
    dx = length / mm
    dy = height / nn
    dt = time / istep_max

    CFL_no = inlet_velocity * dt / dx
    Peclet_no = inlet_velocity * dx / nue

    print(f"\nm, n = {m}, {n}")
    print(f"dx, dy = {dx}, {dy}")
    print(f"dt = {dt}")
    print(f"CFL_no = {CFL_no}")
    print(f"Peclet_no = {Peclet_no}")

    xp[m1] = 0.0  # front of body
    xp[m1 + mm] = length  # back of body
    for i in reversed(range(m1)):
        xp[i] = xp[i + 1] - dx

    for i in range(m1 + 1, m1 + mm):
        xp[i] = xp[i - 1] + dx

    for i in range(m1 + mm + 1, m + 2):
        xp[i] = xp[i - 1] + dx

    yp[n1] = 0.0  # bottom of body
    yp[n1 + nn] = height  # top of body
    for j in reversed(range(n1)):
        yp[j] = yp[j + 1] - dy

    for j in range(n1 + 1, n1 + nn):
        yp[j] = yp[j - 1] + dy

    for j in range(n1 + nn + 1, n + 2):
        yp[j] = yp[j - 1] + dy

    # xp[0] = - 0.5 * dx  # dummy
    # for i in range(1, m + 1):
    #     xp[i] = xp[i - 1] + dx
    # xp[m + 1] = xp[m] + dx  # dummy
    #
    # yp[0] = - 0.5 * dy  # dummy
    # for i in range(1, n + 1):
    #     yp[i] = yp[i - 1] + dy
    # yp[n + 1] = yp[n] + dy  # dummy

    return dx, dy, dt, m, n, istep_max


@jit
def boundary_c(p, u, v, inlet_velocity, outlet_pressure, m, n, xp, yp, height,
               length):
    # inlet (u = inlet_velocity, v = 0, dp/dx = 0 at i = 1)
    for j in range(1, n + 1):
        u[1][j] = inlet_velocity
        v[1][j] = 0.0
        u[0][j] = inlet_velocity  # dummy
        v[0][j] = 0.0  # dummy
        p[0][j] = p[2][j]  # dummy

    # outlet (du/dx = 0, v = 0, p = outlet_pressure at i = m)
    for j in range(1, n + 1):
        u[m + 1][j] = u[m - 1][j]  # dummy
        v[m + 1][j] = 0.0  # dummy
        v[m][j] = 0.0
        p[m][j] = outlet_pressure
        p[m + 1][j] = outlet_pressure  # dummy

    # lower wall (du/dy = 0, v = 0, dp/dy = 0 at j = 1)
    for i in range(m + 2):
        u[i][0] = u[i][2]
        v[i][1] = 0.0
        v[i][0] = 0.0  # dummy
        p[i][0] = p[i][2]  # dummy

    # symmetry (du/dy = 0, v(i,n) = 0, dp/dy = 0, at j = n)
    for i in range(m + 2):
        u[i][n + 1] = u[i][n - 1]
        v[i][n] = 0.0
        v[i][n + 1] = 0.0  # dummy
        p[i][n + 1] = p[i][n - 1]  # dummy

    # inside of body (u = 0, v = 0, p = dummy)
    for i in range(m + 2):
        for j in range(n + 2):
            if (0 <= xp[i] <= length) and (0 <= yp[j] <= height):
                u[i][j] = 0.0
                v[i][j] = 0.0
                p[i][j] = 0.25 * (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] +
                                  p[i][j - 1])

    # # inlet (u(0,j)=inlet_velocity, v=0., p(0,j)=dummy)
    # for j in range(1, n + 1):
    #     u[0][j] = inlet_velocity
    #     v[0][j] = 0.0
    #     p[0][j] = p[1][j]
    #
    # # outlet (du/dx=0., v=0., p=outlet_pressure)
    # for j in range(1, n + 1):
    #     u[m + 1][j] = u[m][j]
    #     v[m + 1][j] = 0.0
    #     p[m + 1][j] = outlet_pressure
    #
    # # lower wall (u_wall(1/lecture10)=0., v(i,0)=0., dp/dy=0.)
    # for i in range(m + lecture10):
    #     u[i][0] = - u[i][1]
    #     v[i][0] = 0.0
    #     p[i][0] = p[i][1]
    #
    # # symmetry  (du/dy=0., v(i,n)=0., v(i,n+1)=dummy, dp/dy=0.)
    # for i in range(m + lecture10):
    #     u[i][n + 1] = u[i][n]  # Staggered grid
    #     v[i][n] = 0.0
    #     v[i][n + 1] = 0.0
    #     p[i][n + 1] = p[i][n]

    return


@jit
def initial_c(p, u, v, inlet_velocity, outlet_pressure, m, n):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            u[i][j] = inlet_velocity
            v[i][j] = 0.0
            p[i][j] = outlet_pressure

    return


@jit
def matrix_c(p, Ap, Ae, Aw, An, As, bb, m, n, xp, yp, height, length):
    # inlet (dp/dx = 0 at i = 1)
    for j in range(1, n + 1):
        Ae[1][j] = Ae[1][j] + Aw[1][j]
        Aw[1][j] = 0.0
        # Ap[1][j] = Ap[1][j] + Aw[1][j]
        # Aw[1][j] = 0.0

    # outlet (p = outlet_pressure at i = m)
    for j in range(1, n + 1):
        bb[m][j] = Ap[m][j] * p[m][j]
        Ae[m][j] = 0.0
        Aw[m][j] = 0.0
        An[m][j] = 0.0
        As[m][j] = 0.0
        # bb[m][j] = bb[m][j] + Ae[m][j] * p[m + 1][j]
        # Ae[m][j] = 0.0

    # lower wall (dp/dy = 0 at j = 1)
    for i in range(1, m + 1):
        An[i][1] = An[i][1] + As[i][1]
        As[i][1] = 0.0
        # Ap[i][1] = Ap[i][1] + As[i][1]
        # As[i][1] = 0.0

    # upper wall (dp/dy = 0 at j = n)
    for i in range(1, m + 1):
        As[i][n] = As[i][n] + An[i][n]
        An[i][n] = 0.0
        # Ap[i][n] = Ap[i][n] + An[i][n]
        # An[i][n] = 0.0

    # wall of body (dp/dx = 0 or dp/dy = 0)
    for i in range(m + 2):
        for j in range(n + 2):
            if xp[i] == 0.0 and (0.0 <= yp[j] <= height):  # front of body
                Ae[i - 1][j] = 0.0

            elif xp[i] == length and (0.0 <= yp[j] <= height):  # back of body
                Aw[i + 1][j] = 0.0

            elif yp[j] == 0.0 and (0.0 <= xp[i] <= length):  # bottom of body
                An[i - 1][j] = 0.0

            elif yp[j] == height and (0.0 <= xp[i] <= length):  # top of body
                As[i + 1][j] = 0.0

    # inside of body (p = dummy)
    for i in range(m + 2):
        for j in range(n + 2):
            if (0 <= xp[i] <= length) and (0 <= yp[j] <= height):
                bb[i][j] = 0.0

    return
