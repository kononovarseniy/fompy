from math import pi

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from fompy.constants import eV, me
from fompy.models import KronigPenneyModel, DiracCombModel, PeriodicPotentialModel
from fompy.units import unit

matplotlib.rc('axes.formatter', useoffset=False)

if __name__ == '__main__':
    # warnings.simplefilter('error')

    def get_band(model: KronigPenneyModel, m):
        ks = np.linspace(0, pi / model.period, 200)

        es = np.linspace(model.u_min-0.1*eV, 0.1*eV, 10000000)
        vs = np.vectorize(model.equation_left_part)(es, m)
        ks2 = np.arccos(np.clip(vs, -1, 1))
        plt.plot(ks2, es/eV)
        plt.show()

        # plt.axvline(model.u_min, color='k', linestyle='--')
        # plt.axhline(-1, color='k', linestyle='--')
        # plt.axhline(1, color='k', linestyle='--')
        # plt.show()

        # es = np.linspace(-0.02*eV, 0, 10000000)
        # ks2 = np.vectorize(model._equation_left)(m, es)
        # plt.plot(es, ks2)
        # plt.ylim(-3/2, 3/2)
        # plt.axhline(-1, color='k', linestyle='--')
        # plt.axhline(1, color='k', linestyle='--')
        # plt.show()

        bands = []
        for band in range(3):
            r = model.estimate_band_range(m, e_coarse=1e-3 * eV, band=band)

            @np.vectorize
            def energy_func(k):
                return model.get_energy(m, k, r, 1e-7 * eV)

            bands.append(energy_func(ks))

        return ks * model.period, bands

    m = 0.49 * me
    U0 = 0.58 * eV
    a = 5 * unit('nm')
    b = 10 * unit('nm')

    kp_model = KronigPenneyModel(a, b, U0)
    #dc_model = DiracCombModel(a + b, a * U0)

    kp_ks, kp_bands = get_band(kp_model, m)

    for i, b in enumerate(kp_bands):
        print(f'Kronig-Penney E{i}: {b[0] / eV}')
        plt.plot(kp_ks, b / eV, label=f'Kronig-Penney #{i}')

    plt.legend()
    plt.minorticks_on()
    plt.xlim(0, pi)
    #plt.axhline(-U0/eV, color='k', linestyle='--')
    plt.grid(True, which='major', axis='both', color='dimgray')
    plt.grid(True, which='minor', axis='both')
    plt.show()
