from math import pi

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from fompy.constants import eV, me
from fompy.models import KronigPenneyModel, DiracCombModel
from fompy.units import unit

matplotlib.rc('axes.formatter', useoffset=False)

if __name__ == '__main__':
    m = 0.49 * me
    U0 = -0.58 * eV  # '-' потенцияальная яма
    a = 1 * unit('nm')
    b = 200 * unit('nm')

    kp_model = KronigPenneyModel(a, b, U0)
    dc_model = DiracCombModel(a + b, a * U0)

    es = np.linspace(-0.001 * eV, 0.001 * eV, 100000)

    ks = kp_model.get_ks(es, m)  # Array of k * (a+b)
    plt.plot(ks, es / eV, label='Kronig-Penney')

    ks = dc_model.get_ks(es, m)
    plt.plot(ks, es / eV, label='Dirac comb')

    plt.axhline(0, color='k', linestyle='--')
    plt.xlim(0, pi)
    plt.xlabel("k*(a+b)")
    plt.ylabel("Energy")
    plt.legend()
    plt.show()
