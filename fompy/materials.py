"""
This module contains classes useful for calculating properties of semiconducting materials.
In addition, it includes instances of those classes for particular materials.

Classes
-------
Semiconductor

DopedSemiconductor : extends Semiconductor

Metal

Objects
-------
Si : Semiconductor
    Silicon at 300 K.
Ge : Semiconductor
    Germanium at 300 K.
GaAs : Semiconductor
    Gallium arsenide at 300 K in Gamma valley.
"""

from functools import partial

from scipy.optimize import bisect

from fompy.constants import me, eV, k
from fompy import phys
from fompy.util.fermi_dirac import fd1


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
    conductivity_type(*, T=None, Ef=None)
        Tell the conductivity type.
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

    def n_concentration(self, Ef=None, T=300):
        """
        Calculate the electron concentration.

        Parameters
        ----------
        Ef=None : float or None
            The Fermi level. If None, `Ef` is found via the `fermi_level` method.
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The electron concentration.
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nc(T) * fd1((Ef - self.Eg) / (k * T))

    def p_concentration(self, Ef=None, T=300):
        """
        Calculate the hole concentration.

        Parameters
        ----------
        Ef=None : float or None
            The Fermi level. If None, `Ef` is found via the `fermi_level` method.
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The hole concentration.
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nv(T) * fd1(-Ef / (k * T))

    def _charge_imbalance(self, Ef, T):
        return self.p_concentration(Ef, T) - self.n_concentration(Ef, T)

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

        Notes
        -----
        All energies are counted from the valence band edge Ev.
        """
        # TODO: Not sure if it works correctly with Nd!=0 and Na!=0
        return bisect(partial(self._charge_imbalance, T=T), 0, self.Eg, xtol=1e-6 * self.Eg)

    def conductivity_type(self, *, T=None, Ef=None):
        """Tell the conductivity type."""
        return 'i'


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

    Methods
    -------
    p_donor_concentration(Ef=None, T=300)
        Calculate the concentration of positive donor ions.
    n_acceptor_concentration(Ef=None, T=300)
        Calculate the concentration of negative acceptor ions.
    conductivity_type(*, T=None, Ef=None)
        Tell the conductivity type (overrides `Semiconductor`).
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
        return self.p_concentration(Ef, T) + self.p_donor_concentration(Ef, T) \
               - self.n_concentration(Ef, T) - self.n_acceptor_concentration(Ef, T)

    def conductivity_type(self, *, T=None, Ef=None):
        """Tell the conductivity type (overrides `Semiconductor`)."""
        if Ef is not None and T is not None:
            raise ValueError('Both T and Ef are specified')
        if T is None:
            T = 300
        if Ef is None:
            Ef = self.fermi_level(T)
        return 'p' if self.p_concentration(Ef, T) > self.n_concentration(Ef, T) else 'n'


class Metal:
    """
    A class to calculate properties of a metal.

    Attributes
    ----------
    work_function : float
        The work function.
    """
    def __init__(self, work_function):
        """
        Construct the necessary attributes for the `Metal` object.

        Parameters
        ----------
        work_function : float
            The work function.
        """
        self.work_function = work_function


# Values at 300K
# http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html
Si = Semiconductor(0.36 * me, 0.81 * me, 1.12 * eV, 4.05 * eV)
"""A `Semiconductor` object for silicon at 300 K."""

# http://www.ioffe.ru/SVA/NSM/Semicond/Ge/index.html
Ge = Semiconductor(0.22 * me, 0.34 * me, 0.661 * eV, 4.0 * eV)
"""A `Semiconductor` object for germanium at 300 K."""

# http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/index.html
GaAs = Semiconductor(0.063 * me, 0.53 * me, 1.424 * eV, 4.07 * eV)  # Gamma-valley
"""A `Semiconductor` object for gallium arsenide at 300 K in Gamma valley."""

# TODO: add more materials
