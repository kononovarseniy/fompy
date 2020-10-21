from math import exp

from fompy.constants import k
from fompy.util.fermi_dirac import fd1  # noqa


def fermi(E, Ef, T):
    return 1 / (1 + exp((E - Ef) / (k * T)))
