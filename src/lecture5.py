import numpy as np
from matplotlib import pyplot as plt

# setup condition
dx = 0.05
x = np.zeros(21)
for i in range(21):
    x[i] = dx * i
f = np.zeros((13, 21))
f[0][10] = 1.  # I.C.
d = 0.7

if __name__ == '__main__':
    for i in range(12):
        for j in range(1, 20):
            f[i + 1][j] = f[i][j] + d * (f[i][j + 1] - 2 * f[i][j] + f[i][j - 1])

    # graph
    fig, ax = plt.subplots()
    p1 = plt.plot(x, f[0], marker="o", clip_on=False, label="t = 0")
    p2 = plt.plot(x, f[3], marker="+", ms=8, clip_on=False, label="t = 0.03")
    p3 = plt.plot(x, f[6], marker="s", clip_on=False, label="t = 0.06")
    p4 = plt.plot(x, f[9], marker="D", clip_on=False, label="t = 0.09")
    p5 = plt.plot(x, f[-1], marker="x", clip_on=False, label="t = 0.12")
    plt.xticks([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"])
    ax.set_xticks(np.linspace(0, 1, 11), minor=True)
    plt.yticks([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"])
    ax.set_yticks(np.linspace(0, 1, 3), minor=True)
    plt.title(f"d = {d}")
    plt.legend(loc='upper right')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('f')
    plt.show()
