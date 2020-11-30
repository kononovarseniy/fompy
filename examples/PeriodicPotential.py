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

    def get_band(model: PeriodicPotentialModel, m):
        ks = np.linspace(0, pi / model.period, 200)

        bands = []
        for band in range(1):
            r = model.estimate_band_range(m, e_coarse=1e-3 * eV, band=band)

            @np.vectorize
            def energy_func(k):
                return model.get_energy(m, k, r, 1e-7 * eV)

            bands.append(energy_func(ks))

        return ks * model.period, bands

    m = 0.49 * me
    U0 = 0.58 * eV
    a = 10 * unit('nm')
    b = 100 * unit('nm')

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
