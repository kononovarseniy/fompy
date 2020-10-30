from math import *

import matplotlib.pyplot as plt
import numpy as np

from fompy.constants import *
from fompy.materials import Si


@np.vectorize
def Efi(mat, T):
    return mat.Eg / 2 + 3 / 4 * k * T * log(mat.mh / mat.me)


if __name__ == '__main__':
    ts = np.linspace(10, 5000, 100)

    print(f'Function argument is -4 at T={Si.Eg / (8 * k)}')
    print(f'Function argument is -1 at T={Si.Eg / (2 * k)}')

    ef = np.vectorize(Si.fermi_level)(ts) / eV
    efi = Efi(Si, ts) / eV

    plt.title('Fermi level')
    plt.plot(ts, efi, label='non-degenerate')
    plt.plot(ts, ef, label='generic')
    plt.xlabel('T')
    plt.ylabel('$E_f$')
    print(np.max(np.abs(efi - ef)))
    plt.legend()
    plt.show()
