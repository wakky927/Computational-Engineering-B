import os

import numpy as np
from matplotlib import pyplot as plt

from func import diff, setup_conditions, fit


if __name__ == '__main__':
    # setup conditions
    print("刻み幅を入力してください\n")
    dt = float(input(">> "))

    n = 100

    t = np.zeros(n + 1)
    f_exact = np.zeros(n + 1)
    f_approx = np.zeros(n + 1)

    setup_conditions.lec3(n, dt, t, f_exact)

    # calculation
    ''' Euler method '''
    f_approx[0], f_approx[1] = f_exact[0], f_exact[1]  # I.C.
    diff.lec3_euler(n, dt, f_approx)
    error = np.abs(f_exact - f_approx)

    # write csv file
    m_e = np.transpose(np.stack([t[:-1], f_exact[:-1], f_approx[:-1], error[:-1]]))
    os.makedirs('../data/lecture3/csv/Euler', exist_ok=True)
    np.savetxt('../data/lecture3/csv/Euler/dt_{}.csv'.format(dt),
                m_e, delimiter=',', fmt='%.16f', header="t,f_exact,f_approx,error")

    # output commandline
    print("\n\nEuler method\n")
    print("time, f_exact, f_approx, error\n")
    for i in range(1, n):
        print("{}, {}, {}, {}".format(t[i], f_exact[i],
                                      f_approx[i], error[i]))

    # graph
    fig, ax = plt.subplots()
    ax.plot(t[:-1], f_exact[:-1], label='exact')
    ax.plot(t[:-1], f_approx[:-1], label='approx_Euler')
    ax.plot(t[:-1], error[:-1], label='error')
    plt.title("Euler method")
    plt.legend(loc='lower right')
    plt.xlim(np.min(t), np.max(t))
    plt.yscale("log")
    # plt.show()
    os.makedirs('../data/lecture3/images/Euler', exist_ok=True)
    fig.savefig('../data/lecture3/images/Euler/dt_{}.png'.format(dt))

    ''' Adams-Bashforth method '''
    f_approx[0], f_approx[1] = f_exact[0], f_exact[1]  # I.C.
    diff.lec3_ab(n, dt, f_approx)
    error = f_exact - f_approx

    # write csv file
    m_ab = np.transpose(np.stack([t[:-1], f_exact[:-1], f_approx[:-1], error[:-1]]))
    os.makedirs('../data/lecture3/csv/Adams-Bashforth', exist_ok=True)
    np.savetxt('../data/lecture3/csv/Adams-Bashforth/dt_{}.csv'.format(dt),
               m_ab, delimiter=',', fmt='%.16f', header="t,f_exact,f_approx,error")

    # output commandline
    print("\n\nAdams-Bashforth method\n")
    print("time, f_exact, f_approx, error\n")
    for i in range(1, n):
        print("{}, {}, {}, {}".format(t[i], f_exact[i],
                                      f_approx[i], error[i]))

    # graph
    fig, ax = plt.subplots()
    ax.plot(t[:-1], f_exact[:-1], label='exact')
    ax.plot(t[:-1], f_approx[:-1], label='approx_Adams-Bashforth')
    ax.plot(t[:-1], error[:-1], label='error')
    plt.title("Adams-Bashforth method")
    plt.legend(loc='lower right')
    plt.xlim(np.min(t), np.max(t))
    plt.yscale("log")
    # plt.show()
    os.makedirs('../data/lecture3/images/Adams-Bashforth', exist_ok=True)
    fig.savefig('../data/lecture3/images/Adams-Bashforth/dt_{}.png'.format(dt))

    ''' Trapezoidal method '''
    f_approx[0], f_approx[1] = f_exact[0], f_exact[1]  # I.C.
    diff.lec3_tra(n, dt, f_approx)
    error = f_exact - f_approx

    # write csv file
    m_tra = np.transpose(np.stack([t[:-1], f_exact[:-1], f_approx[:-1], error[:-1]]))
    os.makedirs('../data/lecture3/csv/Trapezoidal', exist_ok=True)
    np.savetxt('../data/lecture3/csv/Trapezoidal/dt_{}.csv'.format(dt),
               m_tra, delimiter=',', fmt='%.16f', header="t,f_exact,f_approx,error")

    # output commandline
    print("\n\nTrapezoidal method\n")
    print("time, f_exact, f_approx, error\n")
    for i in range(1, n):
        print("{}, {}, {}, {}".format(t[i], f_exact[i],
                                      f_approx[i], error[i]))

    # graph
    fig, ax = plt.subplots()
    ax.plot(t[:-1], f_exact[:-1], label='exact')
    ax.plot(t[:-1], f_approx[:-1], label='approx_Trapezoidal')
    ax.plot(t[:-1], error[:-1], label='error')
    plt.title("Trapezoidal method")
    plt.legend(loc='lower right')
    plt.xlim(np.min(t), np.max(t))
    plt.yscale("log")
    # plt.show()
    os.makedirs('../data/lecture3/images/Trapezoidal', exist_ok=True)
    fig.savefig('../data/lecture3/images/Trapezoidal/dt_{}.png'.format(dt))
