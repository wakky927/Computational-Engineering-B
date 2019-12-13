import os

import numpy as np
from matplotlib import pyplot as plt

from func import diff as df
from func import setup_conditions


if __name__ == '__main__':

    # setup conditions
    print("分割数を決定します.\n最大分割数を入力して下さい.(10の倍数, 0 ~ 20000)")
    print("計算は入力された数の0.1倍刻みで, 計10回行われます.\nex) 100 -> 10, 20, 30, 40, 50, 60, 70, 80, 90, 100")
    N = int(input(">> "))

    if N < 0:
        print("正の整数を入力して下さい.")

    elif N > 20000:
        print("入力された数が大きすぎます.")

    elif N % 10 != 0:
        print("入力された数が10の倍数ではありません.")

    else:
        for n in range(int(0.1 * N), N + 1, int(0.1 * N)):
            print("\n\n分割数: {}".format(n))

            dx = 1 / n
            print("dx: {}\n".format(dx))

            x = np.zeros(n + 1)
            f = np.zeros(n + 1)
            df_approx_f = np.zeros(n + 1)
            df_approx_c = np.zeros(n + 1)
            df_exact = np.zeros(n + 1)

            setup_conditions.lec1(n, dx, x, f, df_exact)

            # 1st derivative by 1st-order forward scheme
            df.first_forward(n, dx, f, df_approx_f)
            abs_f = df_approx_f - df_exact

            # write csv file
            m1 = np.stack([x[1:-1], df_exact[1:-1], df_approx_f[1:-1], abs_f[1:-1]])
            os.makedirs('../data/lecture1/csv/first_forward/{}'.format(N), exist_ok=True)
            np.savetxt('../data/lecture1/csv/first_forward/{}/first_forward_{}.csv'.format(N, n),
            m1, delimiter=',', fmt='%.10f')

            # output commandline
            print("first-order forward scheme\n")
            print("x(i),  df_exact(i),  df_approx(i),  abs(error)")
            for i in range(1, n):
                print("{}, {}, {}, {}".format(x[i], df_exact[i],
                                              df_approx_f[i], abs(df_approx_f[i] - df_exact[i])))

            # each graph
            # fig = plt.figure(2 * n - 1)
            plt.scatter(x[1:-1], df_exact[1:-1], s=20, c='m', marker='.', label='f\'_exact')
            plt.scatter(x[1:-1], df_approx_f[1:-1], s=20, c='y', marker='.', label='f\'_approx')
            plt.scatter(x[1:-1], abs_f[1:-1], s=20, c='c', marker='.', label='f\'_error')
            plt.title("first-order forward scheme divided by {}".format(n))
            plt.legend(loc='right')
            # plt.show()
            os.makedirs('../data/lecture1/images/first_forward/{}'.format(N), exist_ok=True)
            fig.savefig('../data/lecture1/images/first_forward/{}/first_forward_{}.png'.format(N, n))

            # errors graph
            fig = plt.figure(1)
            plt.scatter(x[1:-1], abs_f[1:-1], s=20, marker='.', label='f\'_error_{}'.format(n))
            plt.title("first-order forward scheme")
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
            plt.subplots_adjust(right=0.7)
            plt.yscale("log")
            # plt.show()
            os.makedirs('../data/lecture1/images/errors/first_forward/{}'.format(N), exist_ok=True)
            fig.savefig('../data/lecture1/images/errors/first_forward/{}/first_forward_{}.png'.format(N, n))

            # 1st derivative by 2nd-order central scheme
            df.second_central(n, dx, f, df_approx_c)
            abs_c = df_approx_c - df_exact

            # write csv file
            m2 = np.stack([x[1:-1], df_exact[1:-1], df_approx_c[1:-1], abs_c[1:-1]])
            os.makedirs('../data/lecture1/csv/second_central/{}'.format(N), exist_ok=True)
            np.savetxt('../data/lecture1/csv/second_central/{}/second_central_{}.csv'.format(N, n),
            m2, delimiter=',', fmt='%.10f')

            # output commandline
            print("\nsecond-order central scheme\n")
            print("x(i),  df_exact(i),  df_approx(i),  abs(error)")
            for i in range(1, n):
                print("{}, {}, {}, {}".format(x[i], df_exact[i],
                                              df_approx_c[i], abs(df_approx_c[i] - df_exact[i])))

            # each graph
            fig = plt.figure(2 * n)
            plt.scatter(x[1:-1], df_exact[1:-1], s=20, c='m', marker='.', label='f\'_exact')
            plt.scatter(x[1:-1], df_approx_c[1:-1], s=20, c='y', marker='.', label='f\'_approx')
            plt.scatter(x[1:-1], abs_c[1:-1], s=20, c='c', marker='.', label='f\'_error')
            plt.title("second-order central scheme divided by {}".format(n))
            plt.legend(loc='right')
            # plt.show()
            os.makedirs('../data/lecture1/images/second_central/{}'.format(N), exist_ok=True)
            fig.savefig('../data/lecture1/images/second_central/{}/second_central_{}.png'.format(N, n))

            # errors graph
            fig = plt.figure(2)
            plt.scatter(x[1:-1], abs_c[1:-1], s=20, marker='.', label='f\'_error_{}'.format(n))
            plt.title("second-order central scheme")
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
            plt.subplots_adjust(right=0.7)
            plt.yscale("log")
            # plt.show()
            os.makedirs('../data/lecture1/images/errors/second_central/{}'.format(N), exist_ok=True)
            fig.savefig('../data/lecture1/images/errors/second_central/{}/second_central_{}.png'.format(N, n))
