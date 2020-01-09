import os

import numpy as np
from matplotlib import pyplot as plt

import func.operation
import func.setup_conditions
import func.solve_matrix


if __name__ == '__main__':
    # define parameters
    md = 100  # md > grid size(m)
    m = 11  # grid size (m < md)
    dl = 1.0
    dx = dl / (m - 1)
    c = - 2.0  # source term C
    iter_max = 300
    omega_opt = 2 / (1 + np.sqrt(1 - 0.988**2))

    print(f"m = {m}\nL = {dl}\ndx = {dx}\nc = {c}\n")
    print("omega_opt = {}\n\n".format(omega_opt))

    result = np.zeros((iter_max + 1, 5))

    for o in range(1, iter_max + 1):
        # setup matrix (d^2f/dx^2) = b
        x = np.zeros(m + 1)
        ap = np.zeros(m + 1)
        aw = np.zeros(m + 1)
        ae = np.zeros(m + 1)
        bb = np.zeros(m + 1)

        func.setup_conditions.lec4(m, dl, dx, c, x, ap, aw, ae, bb)

        # B.C. preset
        ''' Delecret B.C. p(0):dummy = no use '''
        bb[1] = 0 * ap[1]  # p = 0 at x = 0
        aw[1] = 0
        ae[1] = 0

        # print("Delecret B.C. at x = 0")

        ''' Neumann B.C. p(m + 1): dummy = no use '''
        ap[m] = ap[m] + ae[m]  # dpdx = 0. at x = L
        ae[m] = 0

        # print("Neumann B.C. at x = L\n")

        # I.C. (include B.C. i = 1, m)
        p = np.zeros(m + 1)
        p_exact = np.zeros(m + 1)

        func.operation.lec4_exact(m, dl, c, x, p_exact)

        # SOR method
        relax_factor = o / 100  # SOR relaxation factor
        result[o][0] = relax_factor
        result[o][1], result[o][2], result[o][3], result[o][4] = func.solve_matrix.SOR3(md, p, ap, ae, aw, bb, m, p_exact, relax_factor)

        # output commandline
        # if result[o][1] < 500:
        #     print(f"omega: {relax_factor}")
        #     print(f"SOR iteration no. {result[o][1]}  -- error = {result[o][2]}, {result[o][3]}, {result[o][4]}\n")

        # # write csv file
        # m_e = np.transpose(np.stack([x[1:], p[1:]]))
        # os.makedirs(f'../data/lecture4/csv/output_SOR', exist_ok=True)
        # np.savetxt(f'../data/lecture4/csv/output_SOR/relax_factor_{relax_factor}.csv',
        #            m_e, delimiter=',', fmt='%.16f', header="x,p")

        # graph1
        # fig, ax = plt.subplots()
        # ax.plot(x[1:], p[1:], label='f(x)')
        # ax.plot(x[1:], p_exact[1:], label='exact')
        # plt.title(f"SOR method (relax_factor = {relax_factor})")
        # plt.legend(loc='upper left')
        # plt.xlim(np.min(x), np.max(x))
        # plt.show()
        # os.makedirs('../data/lecture4/images/output_SOR', exist_ok=True)
        # fig.savefig(f'../data/lecture4/images/output_SOR/relax_factor_{relax_factor}.png')

    # graph2
    fig, ax = plt.subplots()
    ax.scatter(result[:, 0], result[:, 1], marker='.', label='result')
    p = plt.plot([omega_opt, omega_opt], [0, 5000], "orange", linestyle='dashed', label='optimized relax factor')
    plt.legend(loc='lower left')
    plt.xlim(0, 2)
    plt.ylim(0, 5000)
    plt.xlabel('relax factor [-]')
    plt.ylabel('iteration number [times]')
    plt.show()
    # os.makedirs('../data/lecture4/images', exist_ok=True)
    # fig.savefig(f'../data/lecture4/images/relax_factor.png')
