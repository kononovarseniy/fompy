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
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nc(T) * fdk(0.5, (Ef - self.Eg) / (k * T))

    def p_intrinsic(self, Ef=None, T=300):
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nv(T) * fdk(0.5, -Ef / (k * T))

    def _charge_imbalance(self, Ef, T):
        return self.p_intrinsic(Ef, T) - self.n_intrinsic(Ef, T)

    def fermi_level(self, T=300):
        return bisect(partial(self._charge_imbalance, T=T), 0, self.Eg, xtol=1e-6 * self.Eg)


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


# TODO: add classes for other material types, such as metals

# Values at 300K
# http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html
Si = Semiconductor(0.36 * me, 0.81 * me, 1.12 * eV, 4.05 * eV)
# TODO: add more materials
