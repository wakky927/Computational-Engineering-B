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

    print(f"m = {m}\nL = {dl}\ndx = {dx}\nc = {c}\n")

    for omega in range(1, 201):
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

        print("Neumann B.C. at x = L\n")

        # I.C. (include B.C. i = 1, m)
        p = np.zeros(m + 1)
        p_exact = np.zeros(m + 1)

        func.operation.lec4_exact(m, dl, c, x, p_exact)

        # SOR method
        relax_factor = omega / 100  # SOR relaxation factor
        iter_n, error1, error2, error3 = func.solve_matrix.SOR(md, p, ap, ae, aw, bb, m, p_exact, relax_factor)
        print(f"\nSOR iteration no. {iter_n}  -- error = {error1}, {error2}, {error3}\n")

        # output commandline
        print(f"m = {m}\nx    , p")
        for i in range(1, m + 1):
            print(x[i], p[i], i)

        # write csv file
        m_e = np.transpose(np.stack([x[1:], p[1:]]))
        os.makedirs(f'../data/lecture4/csv/output_SOR', exist_ok=True)
        np.savetxt(f'../data/lecture4/csv/output_SOR/relax_factor_{relax_factor}.csv',
                   m_e, delimiter=',', fmt='%.16f', header="x,p")

        # graph
        fig, ax = plt.subplots()
        ax.plot(x[1:], p[1:], label='f(x)')
        ax.plot(x[1:], p_exact[1:], label='exact')
        plt.title(f"SOR method (relax_factor = {relax_factor})")
        plt.legend(loc='upper left')
        plt.xlim(np.min(x), np.max(x))
        # plt.show()
        os.makedirs('../data/lecture4/images/output_SOR', exist_ok=True)
        fig.savefig(f'../data/lecture4/images/output_SOR/relax_factor_{relax_factor}.png')
