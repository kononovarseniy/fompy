from math import exp, pi, sqrt
from fompy.constants import e, k, h_bar


# w depletion area width
def depletion_width(eps, n, d_phi):
    """Ширина зоны между металлом и полупроводником"""
    return sqrt(eps * d_phi / (2 * pi * e * n))


# L_D screening length
def debye_length(eps, n, T):
    """Длина Дебая"""
    return sqrt(eps * k * T / (4 * pi * e ** 2 * n))


def fermi(E, Ef, T):
    return 1 / (1 + exp((E - Ef) / (k * T)))


def conductivity(n, n_mob, p, p_mob):
    return e * (n * n_mob + p *p_mob)


def effective_state_density(m_eff, T):
    """
    Used to calculate Nc and Nv
    
    m_eff -- Effective mass of density of states
    T -- Temperature
    """
    return 2 * (2 * pi * m_eff * k * T / (2 * pi * h_bar) ** 2) ** (3 / 2)

