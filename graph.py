import sys

import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    x = np.array([], dtype=np.float64)
    y = np.array([], dtype=np.float64)
    z = np.array([], dtype=np.float64)

    # filename = '_'.join(next(sys.stdin).split())

    for line in sys.stdin:
        x_value, y_value, z_value = map(float, line.split())
        x = np.append(x, x_value)
        y = np.append(y, y_value)
        z = np.append(z, z_value)

    fig, ax = plt.subplots()

    ax.plot(x, y)
    ax.plot(x, z)
    plt.savefig(f'fig.jpg')
