import matplotlib.pyplot as plt
import numpy as np

from fompy.constants import eV
from fompy.materials import Si
from fompy.models import DopedSemiconductor

if __name__ == '__main__':
    # Хотя материалы 1 и 2 имеют разные концентрации и донорные уровни,
    # разницы между концентрациями доноров и акцепторов одинаковы.
    # Поэтому кривые дисбаланса, а следовательно и уровни ферми, почти совпадают.
    # При совпадении донорных уровней уровни ферми становятся еще ближе.
    # Для вырожденного полупроводника это приближение уже не работает.
    # Это связанно с тем что интеграл Ферми-Дирака не обладает
    # основным функциональным свойством экспоненты exp(a+b) = exp(a)*exp(b).
    # А в выровжденном мы как раз больше не можем приближать интеграл Ферми-Дирака экспонентой.
    mat1 = DopedSemiconductor(Si, 2e17, 0.045 * eV, 1e17, Si.Eg - 0.045 * eV)
    mat2 = DopedSemiconductor(Si, 1e17, 0, 0, Si.Eg)
    mat3 = DopedSemiconductor(Si, 2e17, 0.045 * eV, 0, Si.Eg - 0.045 * eV)
    print('E_f_1 =', mat1.fermi_level() / eV, 'eV')
    print('E_f_2 =', mat2.fermi_level() / eV, 'eV')
    min_e = 0.1 * Si.Eg
    max_e = 0.9 * Si.Eg
    xs = np.linspace(min_e, max_e, 1000)
    func = np.vectorize(DopedSemiconductor._charge_imbalance)
    ys1 = func(mat1, xs, T=300)
    ys2 = func(mat2, xs, T=300)
    ys3 = func(mat3, xs, T=300)

    plt.plot(xs / eV, ys1, 'r', label='1')
    plt.plot(xs / eV, ys2, 'g', label='2')
    plt.plot(xs / eV, ys3, 'b', label='3')

    plt.title('Charge carrier imbalance')
    plt.xlabel('$E_f, eV$')
    plt.ylabel(r'$\Delta N$')
    plt.axvline(x=0)
    plt.axvline(x=Si.Eg / eV)
    plt.grid(b=True, which='both', axis='both')
    plt.legend()
    plt.show()
