from math import exp, pi
from constants import e, k, h_bar
from scipy.optimize import bisect
from fdint import fdk

# w depletion area width
def depletion_width(eps, n, d_phi):
    """Ширина зоны между металлом и полупроводником"""
    return sqrt(eps*d_phi/(2*pi*e*n))

# L_D screening length
def debye_length(eps, n, T):
    """Длина Дебая"""
    return sqrt(eps*k*T/(4*pi*e**2*n))

def fermi(E, Ef, T):
    return 1 / (1 + exp((E - Ef) / (k * T)))

def effective_state_density(m_eff, T):
    """
    Used to calculate Nc and Nv
    
    m_eff -- Effective mass of density of states
    T -- Temperature
    """
    return 2*(2*pi*m_eff*k*T/(2*pi*h_bar)**2)**(3/2)


def n_intrinsic(m_eff, Ec, Ef, T):
    """
    Intrinsic electron concentration
    
    Ef -- The Fermi_level
    Ec -- The conduction band edge
    Ev == 0 -- The valence band edge
    """
    # fdk is Fermi-Dirak integral
    return effective_state_density(m_eff, T) * fdk(0.5, (Ef - Ec) / (k * T))

def p_intrinsic(m_eff, Ef, T):
    """
    Intrinsic hole concentration
    
    Ef -- The Fermi_level
    Ev == 0 -- The valence band edge
    """
    # fdk is Fermi-Dirak integral
    return effective_state_density(m_eff, T) * fdk(0.5, -Ef / (k * T))

def p_donor_concentration(Nd, Ed, Eg, Ef, T):
    return Nd * (1 - fermi(Ed, Ef, T))
    
def n_acceptor_concentration(Na, Ea, Ef, T):
    return Na * fermi(Ea, Ef, T)
    
def fermi_level(me_eff, mh_eff, Eg, Na, Ea, Nd, Ed, T):
    """
    Calculate fermi level (Not shure if it works correctly whith Nd!=0 and Na!=0)
    
    NOTE: all energies are counted from Ev
    """
    def imbalance(Ef):
        return p_intrinsic(mh_eff, Ef, T) + p_donor_concentration(Nd, Ed, Eg, Ef, T) \
            - n_intrinsic(me_eff, Eg, Ef, T) - n_acceptor_concentration(Na, Ea, Ef, T)
    return bisect(imbalance, 0, Eg, xtol=1e-6*Eg)

