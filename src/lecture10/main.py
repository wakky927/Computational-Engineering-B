import numpy as np

import condition
import output
import solver


if __name__ == '__main__':
    md = 202
    nd = 202

    u = np.zeros((md, nd))
    v = np.zeros((md, nd))
    p = np.zeros((md, nd))
    u_old = np.zeros((md, nd))
    v_old = np.zeros((md, nd))

    xp = np.zeros(md)
    yp = np.zeros(nd)

    # setup conditions
    nue, density, length, height, time, inlet_velocity, outlet_pressure\
        = condition.physical_c()
    dx, dy, dt, m, n, istep_max\
        = condition.grid_c(xp, yp, nue, length, height, time, inlet_velocity)
    output.grid(xp, yp, m, n, dt)

    istep = 0
    time = istep * dt

    condition.initial_c(p, u, v, inlet_velocity, outlet_pressure, m, n)
    condition.boundary_c(p, u, v, inlet_velocity, outlet_pressure, m, n,
                         xp, yp, height, length)

    # print initial conditions on commandline
    # output.solution(p, u, v, m, n)

    ''' MAC algorithm start '''
    for istep in range(1, istep_max + 1):
        time = istep * dt
        print(f"--- time_steps = {istep}, -- time = {time}")

        for i in range(m + 2):
            for j in range(n + 2):
                u_old[i][j] = u[i][j]
                v_old[i][j] = v[i][j]

        solver.solve_p(p, u, v, u_old, v_old, nue, density, dx, dy, dt, m, n,
                       xp, yp, height, length)
        solver.solve_u(p, u, density, dx, dt, m, n)
        solver.solve_v(p, v, density, dy, dt, m, n)
        condition.boundary_c(p, u, v, inlet_velocity, outlet_pressure, m, n,
                             xp, yp, height, length)

        # output.solution(p, u, v, m, n)

    output.divergent(p, u, v, dx, dy, m, n, dt)
    ''' MAC algorithm end '''

    # print conditions (recall)
    # condition.physical_c()
    # condition.grid_c(xp, yp, nue, length, height, time, inlet_velocity)

    # print solutions
    # output.solution(p, u, v, m, n)
    output.solution_post(p, u, v, m, n, dt)
    output.paraview(p, xp, yp, m, n, u, v, dt)
