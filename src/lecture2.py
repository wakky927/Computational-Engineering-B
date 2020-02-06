import os

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

from func import operation, setup_conditions, fit


if __name__ == '__main__':
    # setup conditions

    N = 1000
    X = np.zeros(int(N * 0.1) + 1)
    Y = np.zeros(int(N * 0.1) + 1)

    for n in range(int(0.01 * N), N + 1, int(0.01 * N)):
        dx = 1 / n
        X[int(n * 0.1)] = dx

        x = np.zeros(n + 1)
        f = np.zeros(n + 1)

        f1_exact, f2_exact = setup_conditions.lec2(n, dx, x, f)

        f1_approx = operation.lec2_first_func(dx, f)
        f2_approx = operation.lec2_second_func(dx, f)

        # print("\ndx, f1_exact, f1_approx, error")
        # print(dx, f1_exact, f1_approx, abs(f1_approx - f1_exact))

        # print("\ndx, f2_exact, f2_approx, error")
        # print(dx, f2_exact, f2_approx, abs(f2_approx - f2_exact))

        error1 = abs(f1_approx - f1_exact)
        error2 = abs(f2_approx - f2_exact)

        Y[int(n * 0.1)] = error1 - error2

    # graph
    fig, ax = plt.subplots()
    ax.plot(X, Y, lw=0, marker='o', label='Truncation error', clip_on=False)

    ''' 3次曲線近似 '''
    # p = np.poly1d(np.polyfit(X, Y, 3))
    # plt.plot(fit.xp(X), p(fit.xp(X)), label='y = ({})x^3+({})x^lecture10\n      +({})x+({})'.format(p[0], p[1], p[lecture10], p[3]))

    ''' 累乗近似 '''
    opt, cov = optimize.curve_fit(fit.loglog, X[1:], Y[1:])
    ax.plot(fit.xp(X), opt[0]*fit.xp(X)**opt[1], lw=1, label='y = {}x^{}'.format(opt[0], opt[1]))

    plt.title("Error vs. Stepping-width")
    plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', labelspacing=1.25, prop={'size': 8.5, })
    plt.subplots_adjust(bottom=0.3)
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.min(Y), np.max(Y))
    # plt.show()
    os.makedirs('../data/lecture2/images'.format(N), exist_ok=True)
    fig.savefig('../data/lecture2/images/{}.png'.format(N))

    # write csv file
    # m2 = np.transpose(np.stack([X[1:], Y[1:]]))
    # os.makedirs('../data/lecture2/csv', exist_ok=True)
    # np.savetxt('../data/lecture2/csv/{}.csv'.format(N),
    #             m2, delimiter=',', fmt='%.16f', header="dx,error")

