"""
This module contains classes useful for calculating properties of semiconducting materials.
In addition, it includes instances of those classes for particular materials.

Classes
-------
Semiconductor

DopedSemiconductor : extends Semiconductor

Objects
-------
Si : Semiconductor
    Silicon at 300 K.
"""

from functools import partial

from fdint import fdk
from scipy.optimize import bisect

from fompy.constants import me, eV, k
from fompy import phys


class Semiconductor:
    """
    A class to calculate properties of an intrinsic (pure) semiconductor.

    Attributes
    ----------
    me : float
        The effective mass of an electron.
    mh : float
        The effective mass of a hole.
    Eg : float
        The energy gap.
    chi : float
        The electron affinity.

    Methods
    -------
    Nc(T=300)
        Calculate the effective density of states for electrons in the conduction band.
    Nv(T=300)
        Calculate the effective density of states for holes in the valence band.
    n_intrinsic(Ef=None, T=300)
        Calculate the intrinsic electron concentration.
    p_intrinsic(Ef=None, T=300)
        Calculate the intrinsic hole concentration.
    fermi_level(T=300):
        Determine the Fermi level from the condition of electroneutrality.
    """

    def __init__(self, me_eff, mh_eff, Eg, chi):
        """
        Construct the necessary attributes for the `Semiconductor` object.

        Parameters
        ----------
        me_eff : float
            The effective mass of an electron.
        mh_eff : float
            The effective mass of a hole.
        Eg : float
            The energy gap.
        chi : float
            The electron affinity.
        """
        self.me = me_eff
        self.mh = mh_eff
        self.Eg = Eg
        self.chi = chi

    def Nc(self, T=300):
        """
        Calculate the effective density of states for electrons in the conduction band.

        Parameters
        ----------
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The effective density of states for electrons in the conduction band.
        """
        return phys.effective_state_density(self.me, T)

    def Nv(self, T=300):
        """
        Calculate the effective density of states for holes in the valence band.

        Parameters
        ----------
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The effective density of states for holes in the valence band.
        """
        return phys.effective_state_density(self.mh, T)

    def n_intrinsic(self, Ef=None, T=300):
        """
        Calculate the intrinsic electron concentration.

        Parameters
        ----------
        Ef=None : float or None
            The Fermi level. If None, `Ef` is found via the `fermi_level` method.
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The intrinsic electron concentration.
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nc(T) * fdk(0.5, (Ef - self.Eg) / (k * T))

    def p_intrinsic(self, Ef=None, T=300):
        """
        Calculate the intrinsic hole concentration.

        Parameters
        ----------
        Ef=None : float or None
            The Fermi level. If None, `Ef` is found via the `fermi_level` method.
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The intrinsic hole concentration.
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nv(T) * fdk(0.5, -Ef / (k * T))

    def _charge_imbalance(self, Ef, T):
        return self.p_intrinsic(Ef, T) - self.n_intrinsic(Ef, T)

    def fermi_level(self, T=300):
        """
        Determine the Fermi level from the condition of electroneutrality.

        Parameters
        ----------
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The Fermi level.
        """
        return bisect(partial(self._charge_imbalance, T=T), 0, self.Eg, xtol=1e-6 * self.Eg)


class DopedSemiconductor(Semiconductor):
    """
    A class to calculate properties of a doped semiconductor (extends `Semiconductor`).

    Attributes
    ----------
    Na : float
        The acceptor concentration.
    Ea : float
        The acceptor level.
    Nd : float
        The donor concentration.
    Ed : float
        The donor level.
    """

    def __init__(self, mat, Na, Ea, Nd, Ed):
        """
        Construct the necessary attributes for the `DopedSemiconductor` object.

        Parameters
        ----------
        mat : Semiconductor
            The intrinsic (pure) semiconductor.
        Na : float
            The acceptor concentration.
        Ea : float
            The acceptor level.
        Nd : float
            The donor concentration.
        Ed : float
            The donor level.
        """
        super(DopedSemiconductor, self).__init__(mat.me, mat.mh, mat.Eg, mat.chi)
        self.Na = Na
        self.Ea = Ea
        self.Nd = Nd
        self.Ed = Ed

    def p_donor_concentration(self, Ef=None, T=300):
        """
        Calculate the concentration of positive donor ions.

        Parameters
        ----------
        Ef=None : float or None
            The Fermi level. If None, `Ef` is found via the `fermi_level` method.
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The concentration of positive donor ions.
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nd * (1 - phys.fermi(self.Ed, Ef, T))

    def n_acceptor_concentration(self, Ef=None, T=300):
        """
        Calculate the concentration of negative acceptor ions.

        Parameters
        ----------
        Ef=None : float or None
            The Fermi level. If None, `Ef` is found via the `fermi_level` method.
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The concentration of negative acceptor ions.
        """
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
"""A `Semiconductor` object for silicon at 300 K."""
# TODO: add more materials
