import matplotlib.pyplot as plt
import numpy as np

from fompy.constants import eV
from fompy.materials import Si
from fompy.models import DopedSemiconductor

if __name__ == '__main__':
    mat1 = DopedSemiconductor(Si, 1e18, 0.045 * eV, 1e20, Si.Eg - 0.045 * eV)
    mat2 = DopedSemiconductor(Si, 0, 0.045 * eV, 1e15, Si.Eg - 0.045 * eV)
    mat3 = DopedSemiconductor(Si, 1e20, 0.4 * eV, 1e20, Si.Eg - 0.4 * eV)
    mat4 = DopedSemiconductor(Si, 2e20, 0.045 * eV, 0, Si.Eg)
    min_e = -0.1 * Si.Eg
    max_e = 1.1 * Si.Eg
    xs = np.linspace(min_e, max_e, 1000)
    func = np.vectorize(DopedSemiconductor._charge_imbalance)
    ys1 = func(mat1, xs, T=300)
    ys2 = func(mat2, xs, T=300)
    ys3 = func(mat3, xs, T=300)
    ys4 = func(mat4, xs, T=300)

    plt.plot(xs / eV, ys1, 'r', label='1')
    plt.plot(xs / eV, ys2, 'g', label='2')
    plt.plot(xs / eV, ys3, 'b', label='3')
    plt.plot(xs / eV, ys4, 'y', label='4')

    plt.title('Charge carrier imbalance')
    plt.xlabel('$E_f, eV$')
    plt.ylabel(r'$\Delta N$')
    plt.axvline(x=0)
    plt.axvline(x=Si.Eg / eV)
    plt.grid(b=True, which='both', axis='both')
    plt.legend()
    plt.show()
