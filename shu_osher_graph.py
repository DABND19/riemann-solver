import sys
from typing import TextIO, Tuple

import matplotlib.pyplot as plt
import numpy as np


def parse_file(file_: TextIO) -> Tuple[
    np.ndarray, 
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray
]:
    x = np.array([])
    p = np.array([])
    rho = np.array([])
    u = np.array([])
    Mach = np.array([])
    for line in file_:
        x_, p_, rho_, u_, Mach_ = map(float, line.split())
        x = np.append(x, x_)
        p = np.append(p, p_)
        rho = np.append(rho, rho_)
        u = np.append(u, u_)
        Mach = np.append(Mach, Mach_)
    return x, p, rho, u, Mach


if __name__ == '__main__':
    fig, ((ax1, ax2), 
          (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, 
                                     figsize=(10, 8), dpi=200)

    with open('shu_osher_godunov_10000_cells.txt') as f:
        x, p, rho, u, Mach = parse_file(f)
        legend = 'Exact'

        ax1.plot(x, p, '--', label=legend)
        ax2.plot(x, rho, '--', label=legend)
        ax3.plot(x, u, '--', label=legend)
        ax4.plot(x, Mach, '--', label=legend)

    legends = {
        'shu_osher_godunov_500_cells.txt': 'Godunov',
        'shu_osher_hllc_500_cells.txt': 'HLLC',
        'shu_osher_hll_500_cells.txt': 'HLL'
    }
    for filename, legend in reversed(legends.items()):
        with open(filename) as f:
            x, p, rho, u, Mach = parse_file(f)

        ax1.plot(x, p, '--', label=legend)
        ax2.plot(x, rho, '--', label=legend)
        ax3.plot(x, u, '--', label=legend)
        ax4.plot(x, Mach, '--', label=legend)

    ax1.set_xlabel('x')
    ax1.set_ylabel('Pressure')
    ax1.grid(True)
    ax1.legend()
    ax2.set_xlabel('x')
    ax2.set_ylabel('Density')
    ax2.grid(True)
    ax2.legend()
    ax3.set_xlabel('x')
    ax3.set_ylabel('Velocity')
    ax3.grid(True)
    ax3.legend()
    ax4.set_xlabel('x')
    ax4.set_ylabel('Mach number')
    ax4.grid(True)
    ax4.legend()

    plt.savefig('pictures/comparison.png')

