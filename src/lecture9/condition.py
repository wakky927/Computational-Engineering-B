def physical_c():
    nue = 1.0e-4  # [m^lecture10/s] : water
    density = 1000.0  # [kg/m^3] : water

    length = 1.0  # [m]
    height = 0.1  # [m]
    time = 100.0  # [s]
    inlet_velocity = 0.1  # [m/s]
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
    m = 10
    n = 10
    istep_max = 1000
    dx = length / m
    dy = height / n
    dt = time / istep_max

    CFL_no = inlet_velocity * dt / dx
    Peclet_no = inlet_velocity * dx / nue

    print(f"\nm, n = {m}, {n}")
    print(f"dx, dy = {dx}, {dy}")
    print(f"dt = {dt}")
    print(f"CFL_no = {CFL_no}")
    print(f"Peclet_no = {Peclet_no}")

    xp[1] = 0.0
    for i in range(2, m + 1):
        xp[i] = xp[i - 1] + dx
    xp[0] = xp[1] - dx  # dummy
    xp[m + 1] = xp[m] + dx  # dummy

    yp[1] = 0.0
    for i in range(2, n + 1):
        yp[i] = yp[i - 1] + dy
    yp[0] = yp[1] - dy  # dummy
    yp[n + 1] = yp[n] + dy  # dummy

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


def boundary_c(p, u, v, inlet_velocity, outlet_pressure, m, n):
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

    # lower wall (u = 0, v = 0, dp/dy = 0 at j = 1)
    for i in range(m + 2):
        u[i][1] = 0.0
        v[i][1] = 0.0
        u[i][0] = 0.0  # dummy
        v[i][0] = 0.0  # dummy
        p[i][0] = p[i][2]  # dummy

    # symmetry (du/dy = 0, v(i,n) = 0, dp/dy = 0, at j = n)
    for i in range(m + 2):
        u[i][n + 1] = u[i][n - 1]
        v[i][n] = 0.0
        v[i][n + 1] = 0.0  # dummy
        p[i][n + 1] = p[i][n - 1]  # dummy

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


def initial_c(p, u, v, inlet_velocity, outlet_pressure, m, n):
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            u[i][j] = inlet_velocity
            v[i][j] = 0.0
            p[i][j] = outlet_pressure

    return


def matrix_c(p, Ap, Ae, Aw, An, As, bb, m, n):
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

    # symmetry (dp/dy = 0 at j = n)
    for i in range(1, m + 1):
        As[i][n] = As[i][n] + An[i][n]
        An[i][n] = 0.0
        # Ap[i][n] = Ap[i][n] + An[i][n]
        # An[i][n] = 0.0

    return


