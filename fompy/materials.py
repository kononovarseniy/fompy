from functools import partial

from fdint import fdk
from scipy.optimize import bisect

from fompy.constants import me, eV, k
from fompy import phys


class Semiconductor:
    def __init__(self, me_eff, mh_eff, Eg, chi):
        """
        _me -- effective mass of electron
        _mh -- effective mass of hole
        _gap -- energy gap
        _chi -- electron affinity
        """
        self.me = me_eff
        self.mh = mh_eff
        self.Eg = Eg
        self.chi = chi

    def Nc(self, T=300):
        return phys.effective_state_density(self.me, T)

    def Nv(self, T=300):
        return phys.effective_state_density(self.mh, T)

    def n_intrinsic(self, Ef=None, T=300):
        """
        Intrinsic electron concentration
        
        Ef -- The Fermi_level
        Ev == 0 -- The valence band edge
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nc(T) * fdk(0.5, (Ef - self.Eg) / (k * T))

    def p_intrinsic(self, Ef=None, T=300):
        """
        Intrinsic hole concentration
        
        Ef -- The Fermi_level
        Ev == 0 -- The valence band edge
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nv(T) * fdk(0.5, -Ef / (k * T))

    def _charge_imbalance(self, Ef, T):
        return self.p_intrinsic(Ef, T) - self.n_intrinsic(Ef, T)

    def fermi_level(self, T=300):
        """
        Calculate fermi level (Not shure if it works correctly whith Nd!=0 and Na!=0)
        
        NOTE: all energies are counted from Ev
        """
        return bisect(partial(self._charge_imbalance, T=T), 0, self.Eg, xtol=1e-6 * self.Eg)

    def conductivity_type(self, *, T=None, Ef=None):
        return 'i'


class DopedSemiconductor(Semiconductor):
    def __init__(self, mat, Na, Ea, Nd, Ed):
        super(DopedSemiconductor, self).__init__(mat.me, mat.mh, mat.Eg, mat.chi)
        self.Na = Na
        self.Ea = Ea
        self.Nd = Nd
        self.Ed = Ed

    def p_donor_concentration(self, Ef=None, T=300):
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nd * (1 - phys.fermi(self.Ed, Ef, T))

    def n_acceptor_concentration(self, Ef=None, T=300):
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Na * phys.fermi(self.Ea, Ef, T)

    def _charge_imbalance(self, Ef, T):
        return self.p_intrinsic(Ef, T) + self.p_donor_concentration(Ef, T) \
               - self.n_intrinsic(Ef, T) - self.n_acceptor_concentration(Ef, T)

    def conductivity_type(self, *, T=None, Ef=None):
        if Ef is not None and T is not None:
            raise ValueError('Both T and Ef are specified')
        if T is None:
            T = 300
        if Ef is None:
            Ef = self.fermi_level(T)
        return 'p' if self.p_intrinsic(Ef, T) > self.n_intrinsic(Ef, T) else 'n'


class Metal:
    def __init__(self, work_function):
        self.work_function = work_function


# Values at 300K
# http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html
Si = Semiconductor(0.36 * me, 0.81 * me, 1.12 * eV, 4.05 * eV)

# http://www.ioffe.ru/SVA/NSM/Semicond/Ge/index.html
Ge = Semiconductor(0.22 * me, 0.34 * me, 0.661 * eV, 4.0 * eV)

# http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/index.html
GaAs = Semiconductor(0.063 * me, 0.53 * me, 1.424 * eV, 4.07 * eV)  # Gamma-valley

# TODO: add more materials
