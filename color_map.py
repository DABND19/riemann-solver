import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


if __name__ == '__main__':
    pictures_dir = Path('pictures/')

    try:
        _, filename = sys.argv
    except ValueError:
        print('You must specify filename',
              file=sys.stderr)
        exit(1)

    data = json.load(sys.stdin)
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
                                   figsize=(10, 4), dpi=200)

    t = np.array(list(map(lambda x: x['t'], data)), np.float64)
    x = np.array(data[0]['x'])

    density = np.array(
        list(map(lambda x: np.array(x['density'], np.float64), data)))
    pressure = np.array(
        list(map(lambda x: np.array(x['pressure'], np.float64), data)))

    density_map = ax1.imshow(
        density[::-1], interpolation='bilinear', norm=LogNorm(vmin=1e-11, vmax=1.5),
        extent=[min(x), max(x), min(t), max(t)], aspect='auto')
    plt.colorbar(density_map, ax=ax1)
    ax1.set_xlabel('Radius')
    ax1.set_ylabel('Density')

    pressure_map = ax2.imshow(
        pressure[::-1], interpolation='bilinear', norm=LogNorm(vmin=1e-11, vmax=1.5),
        extent=[min(x), max(x), min(t), max(t)], aspect='auto')
    plt.colorbar(pressure_map, ax=ax2)
    ax2.set_xlabel('Radius')
    ax2.set_ylabel('Pressure')

    plt.savefig(pictures_dir.joinpath(filename))
