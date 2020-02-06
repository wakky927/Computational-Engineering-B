import os

import numpy as np
from matplotlib import pyplot as plt


if __name__ == '__main__':
    e10 = np.loadtxt('../../data/lecture9/e10.csv', delimiter=',')
    e20 = np.loadtxt('../../data/lecture9/e20.csv', delimiter=',')
    y10 = np.zeros(10)
    y20 = np.zeros(20)

    for i in range(10):
        y10[i] = i * 0.01

    for j in range(20):
        y20[j] = j * 0.005

    fig1, ax1 = plt.subplots(figsize=(8, 6))
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 10
    e1 = ax1.plot(y10[1:], e10, label='10', color='blue')
    e2 = ax1.plot(y20[1:], e20, label='20', color='orange')
    plt.title(f"error (n, m = 10 and 20)")
    ax1.legend(loc='lower right')
    ax1.set_xlabel('y [m]')
    ax1.set_ylabel('error [%]')
    plt.show()
    fig1.savefig(f'../../data/lecture9/1/error_10_20.png')

    print("fin")
