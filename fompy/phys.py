from enum import Enum
from math import exp, pi, sqrt

from fompy.constants import e, k, h_bar


def fermi(E, Ef, T):
    return 1 / (1 + exp((E - Ef) / (k * T)))


def conductivity(n, n_mob, p, p_mob):
    return e * (n * n_mob + p * p_mob)


def effective_state_density(m_eff, T):
    """
    Used to calculate Nc and Nv
    
    m_eff -- Effective mass of density of states
    T -- Temperature
    """
    return 2 * (2 * pi * m_eff * k * T / (2 * pi * h_bar) ** 2) ** (3 / 2)


# w depletion area width
def depletion_width(eps, n, d_phi):
    """Ширина зоны между металлом и полупроводником"""
    return sqrt(eps * d_phi / (2 * pi * e * n))


# L_D screening length
def debye_length(eps, n, T):
    """Длина Дебая"""
    return sqrt(eps * k * T / (4 * pi * e ** 2 * n))


class ContactType(Enum):
    AUGMENTATION = 0,
    DEPLETION = 1,
    INVERSION = 2


class MetalSemiconductorContact:
    def __init__(self, metal, sc):
        self.metal = metal
        self.sc = sc

    def delta_phi(self, T=300):
        """Returns difference of exit potential of semiconductor and metal (in units of voltage)"""
        return (self.sc.Eg - self.sc.fermi_level(T) + self.sc.chi - self.metal.work_function) / e

    def contact_type(self, T=300):
        ct = self.sc.conductivity_type(T=T)
        Ef = self.sc.fermi_level(T)
        dEf = self.delta_phi(T) * e
        if ct == 'i':
            raise NotImplemented('I do not know which cases are possible')
        else:
            if ct != self.sc.conductivity_type(Ef=Ef + dEf):
                return ContactType.INVERSION
            else:
                if ct == 'p':
                    return ContactType.AUGMENTATION if dEf < 0 else ContactType.DEPLETION
                if ct == 'n':
                    return ContactType.DEPLETION if dEf < 0 else ContactType.AUGMENTATION
