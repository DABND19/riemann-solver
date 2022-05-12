import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    pictures_dir = Path('pictures/')

    try:
        _, filename = sys.argv
    except ValueError:
        print('You must specify filename', 
              file=sys.stderr)
        exit(1)

    radius = np.array([], dtype=np.float64)
    pressure = np.array([], dtype=np.float64)
    density = np.array([], dtype=np.float64)
    velocity = np.array([], dtype=np.float64)
    mach_number = np.array([], dtype=np.float64)

    for line in sys.stdin:
        r, p, rho, u, m = map(float, line.split())
        radius = np.append(radius, r)
        pressure = np.append(pressure, p)
        density = np.append(density, rho)
        velocity = np.append(velocity, u)
        mach_number = np.append(mach_number, m)

    fig, ((ax1, ax2), 
          (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, 
                                     figsize=(10, 8), dpi=200)

    ax1.plot(radius, pressure)
    ax1.set_yscale('log')
    ax1.set_xlabel('Radius')
    ax1.set_ylabel('Pressure')

    ax2.plot(radius, density)
    ax2.set_yscale('log')
    ax2.set_xlabel('Radius')
    ax2.set_ylabel('Density')

    ax3.plot(radius, velocity)
    ax3.set_yscale('log')
    ax3.set_xlabel('Radius')
    ax3.set_ylabel('Velocity')

    ax4.plot(radius, mach_number)
    ax4.set_yscale('log')
    ax4.set_xlabel('Radius')
    ax4.set_ylabel('Mach number')

    plt.savefig(pictures_dir.joinpath(filename))
