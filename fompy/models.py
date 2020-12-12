"""
This module contains classes useful for calculating properties of different materials.

Notes
-----
All energies are counted from the valence band edge E<sub>v</sub> &equiv; 0.
"""
import cmath
import math
from abc import ABC, abstractmethod
from enum import Enum
from functools import partial
from math import pi, sqrt, exp, cos

import numpy as np
from scipy.optimize import bisect

from fompy.constants import e, k, h_bar, eV
from fompy.functions import fermi, fd1
from fompy.util.zeros import find_nth_function_zero


def conductivity(n, mobility):
    r"""
    Calculate the conductivity of a material.
    If there are several types of charge carriers, you should sum up the conductivities calculated for these types.

    .. math::
        \sigma = e n \mu

    Parameters
    ----------
    n : float
        The carrier concentration [cm<sup>&minus;3</sup>].
    mobility : float
        The carrier mobility [cm<sup>2</sup> statV<sup>&minus;1</sup> s<sup>&minus;1</sup>].

    Returns
    -------
    float
        The conductivity [s<sup>&minus;1</sup>].
    """
    return e * n * mobility


def concentration(resistivity, mobility):
    r"""
    Calculate the concentration of charge carriers.

    .. math::
        n = \frac{1}{\rho \mu e}

    Parameters
    ----------
    resistivity : float
        The resistivity of the material [s].
    mobility : float
        The carrier mobility [cm<sup>2</sup> statV<sup>&minus;1</sup> s<sup>&minus;1</sup>].

    Returns
    -------
    float
        The carrier concentration [cm<sup>&minus;3</sup>].
    """
    return 1 / (resistivity * mobility * e)


def depletion_width(eps, n, d_phi):
    r"""
    Calculate the width of the depletion region, using the approximation of full depletion.

    .. math::
        w = \sqrt{ \frac{ \epsilon \Delta \phi }{ 2 \pi e n } }

    Parameters
    ----------
    eps : float
        The dielectric constant [1].
    n : float
        The carrier concentration [cm<sup>&minus;3</sup>].
    d_phi : float
        The difference of potentials [statV].
    Returns
    -------
    float
        The width of the depletion region [cm].
    """
    return sqrt(eps * d_phi / (2 * pi * e * n))


def debye_length(eps, n, T):
    r"""
    Calculate the Debye length (screening length).

    .. math::
        \lambda_D = \sqrt{ \frac{ \epsilon k T }{ 4 \pi e^2 n } }

    Parameters
    ----------
    eps : float
        The dielectric constant [1].
    n : float
        The carrier concentration [cm<sup>&minus;3</sup>].
    T : float
        The temperature [K].

    Returns
    -------
    float
        The Debye length [cm].
    """
    return sqrt(eps * k * T / (4 * pi * e ** 2 * n))


def hydrogen_like_energy(eps, m):
    r"""
    Calculate the energy of the electron in a hydrogen-like atom.

    .. math::
        E = \frac{1}{\epsilon^2} \cdot \frac{m_{eff} e^4}{2 {\hbar}^2}

    Parameters
    ----------
    eps : float
        The dielectric constant [1].
    m : float
        The reduced (effective) mass of the electron [g].

    Returns
    -------
    float
        The energy of the electron [erg].
    """
    return e ** 4 * m / (2 * h_bar ** 2) / eps ** 2


def hydrogen_like_radius(eps, m):
    r"""
    Calculate the radius of a hydrogen-like atom.

    ..math::
        r = \frac{\epsilon {\hbar}^2}{m_{eff} e^2}

    Parameters
    ----------
    eps : float
        The dielectric constant [1].
    m : float
        The reduced (effective) mass of the electron [g].

    Returns
    -------
    float
        The radius of the atom [cm].
    """
    return eps * h_bar ** 2 / (m * e ** 2)


class CrystalLattice:
    """
    A class to calculate properties of a crystal lattice.
    """

    def __init__(self, a, m, r, N):
        """
        Construct the necessary attributes for the `CrystalLattice` object.

        Parameters
        ----------
        a : float
            The lattice parameter [cm].
        m : float
            The mass of an atom (the average mass in case of there being several types of atoms) [g].
        r : float
            The maximum radius of non-intersecting spheres around atoms [cm].
        N : float
            The number of atoms in a cell [1].
        """
        self._a = a
        self._m = m
        self._r = r
        self._N = N

    @property
    def a(self):
        """Get the lattice parameter [cm]."""
        return self._a

    @property
    def m(self):
        """Get the mass of an atom (the average mass in case of there being several types of atoms) [g]."""
        return self._m

    @property
    def r(self):
        """Get the maximum radius of non-intersecting spheres around atoms [r]."""
        return self._r

    @property
    def N(self):
        """Get the number of atoms in a cell [1]."""
        return self._N

    @property
    def packing_density(self):
        r"""
        Get the packing density [1].

        .. math::
            \eta = \frac{N}{a^3} \frac{4}{3} \pi r^3
        """
        return self._N * 4 / 3 * pi * self._r ** 3 / self._a ** 3

    @property
    def concentration(self):
        r"""
        Get the concentration [cm<sup>&minus;3</sup>].

        .. math::
            n = \frac{N}{a^3}
        """
        return self._N / self._a ** 3

    @property
    def density(self):
        r"""
        Get the density [g cm<sup>&minus;3</sup>].

        .. math::
            \rho = n m = \frac{N}{a^3} m
        """
        return self.concentration * self._m


class PrimitiveCubicLattice(CrystalLattice):
    """
    A class to calculate properties of a primitive cubic lattice (extends `CrystalLattice`).
    """

    def __init__(self, a, m):
        r"""
        Construct the necessary parameters for the `PrimitiveCubicLattice` object.
        Calls the base class constructor, replacing the following parameters:

        .. math::
            r = \frac{a}{2}
        .. math::
            N = 1

        Parameters
        ----------
        a : float
            The lattice parameter [cm].
        m : float
            The mass of an atom (the average mass in case of there being several types of atoms) [g].
        """
        super().__init__(a, m, a / 2, 1)


class FaceCenteredCubicLattice(CrystalLattice):
    """
    A class to calculate properties of a face-centered cubic lattice (extends `CrystalLattice`).
    """

    def __init__(self, a, m):
        r"""
        Construct the necessary parameters for the `FaceCenteredCubicLattice` object.
        Calls the base class constructor, replacing the following parameters:

        .. math::
            r = \frac{a \sqrt{2}}{4}
        .. math::
            N = 4

        Parameters
        ----------
        a : float
            The lattice parameter [cm].
        m : float
            The mass of an atom (the average mass in case of there being several types of atoms) [g].
        """
        super().__init__(a, m, a * sqrt(2) / 4, 4)


class BodyCenteredCubicLattice(CrystalLattice):
    """
    A class to calculate properties of a body-centered cubic lattice (extends `CrystalLattice`).
    """

    def __init__(self, a, m):
        r"""
        Construct the necessary parameters for the `BodyCenteredCubicLattice` object.
        Calls the base class constructor, replacing the following parameters:

        .. math::
            r = \frac{a \sqrt{3}}{4}
        .. math::
            N = 2

        Parameters
        ----------
        a : float
            The lattice parameter [cm].
        m : float
            The mass of an atom (the average mass in case of there being several types of atoms) [g].
        """
        super().__init__(a, m, a * sqrt(3) / 4, 2)


class DiamondLikeLattice(CrystalLattice):
    """
    A class to calculate properties of a diamond-like lattice.
    """

    def __init__(self, a, m):
        r"""
        Construct the necessary parameters for the `DiamondLikeLattice` object.
        Calls the base class constructor, replacing the following parameters:

        .. math::
            r = \frac{a \sqrt{3}}{8}
        .. math::
            N = 8

        Parameters
        ----------
        a : float
            The lattice parameter [cm].
        m : float
            The mass of an atom (the average mass in case of there being several types of atoms) [g].
        """
        super().__init__(a, m, a * sqrt(3) / 8, 8)


class Semiconductor:
    """
    A class to calculate properties of an intrinsic (pure) semiconductor.

    Attributes
    ----------
    me : float
        The effective mass of an electron [g].
    mh : float
        The effective mass of a hole [g].
    Eg : float
        The energy gap [erg].
    chi : float or None
        The electron affinity [erg].
    eps : float or None
        The dielectric constant [1].
    lattice : CrystalLattice or None
        The crystal lattice.
    """

    def __init__(self, me_eff, mh_eff, Eg, chi=None, eps=None, lattice=None):
        """
        Construct the necessary attributes for the `Semiconductor` object.

        Parameters
        ----------
        me_eff : float
            The effective mass of an electron [g].
        mh_eff : float
            The effective mass of a hole [g].
        Eg : float
            The energy gap [erg].
        chi : float or None
            The electron affinity [erg].
        eps : float or None
            The dielectric constant [1].
        lattice : CrystalLattice or None
            The crystal lattice.
        """
        self.lattice = lattice
        self.me = me_eff
        self.mh = mh_eff
        self.Eg = Eg
        self.chi = chi
        self.eps = eps
        self.donor_energy = dict()
        self.acceptor_energy = dict()

    @staticmethod
    def effective_state_density(m_eff, T):
        r"""
        Calculate the effective density of states.

        .. math::
            N = 2 \left( \frac{2 \pi m_{eff} k T }{ (2 \pi \hbar)^2 } \right)^{3/2}

        Parameters
        ----------
        m_eff : float
            The effective mass [g].
        T : float
            The temperature [K].

        Returns
        -------
        float
            The effective density of states [cm<sup>&minus;3</sup>].
        """
        return 2 * (2 * pi * m_eff * k * T / (2 * pi * h_bar) ** 2) ** (3 / 2)

    def Nc(self, T=300):
        r"""
        Calculate the effective density of states for electrons in the conduction band.

        .. math::
            N_c = 2 \left( \frac{2 \pi m_e k T }{ (2 \pi \hbar)^2 } \right)^{3/2}

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The effective density of states for electrons in the conduction band [cm<sup>&minus;3</sup>].
        """
        return Semiconductor.effective_state_density(self.me, T)

    def Nv(self, T=300):
        r"""
        Calculate the effective density of states for holes in the valence band.

        .. math::
            N_c = 2 \left( \frac{2 \pi m_h k T }{ (2 \pi \hbar)^2 } \right)^{3/2}

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The effective density of states for holes in the valence band [cm<sup>&minus;3</sup>].
        """
        return Semiconductor.effective_state_density(self.mh, T)

    def i_concentration(self, T=300):
        """
        Calculate the intrinsic concentration of electrons (independent of doping).

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The intrinsic concentration of electrons [cm<sup>&minus;3</sup>].
        """
        return self.n_concentration(self.intrinsic_fermi_level(), T)

    def n_concentration(self, Ef=None, T=300):
        r"""
        Calculate the electron concentration.

        .. math::
            n_e = N_c(T) \Phi_{1/2}\left( \frac{ E_f - E_g }{ k T } \right)

        Parameters
        ----------
        Ef : float or None
            The Fermi level [erg]. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
            The temperature [K].

        Returns
        -------
        float
            The electron concentration [cm<sup>&minus;3</sup>].
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nc(T) * fd1((Ef - self.Eg) / (k * T))

    def p_concentration(self, Ef=None, T=300):
        r"""
        Calculate the hole concentration.

        .. math::
            n_h = N_v(T) \Phi_{1/2}\left( \frac{ - E_f }{ k T } \right)

        Parameters
        ----------
        Ef : float or None
            The Fermi level [erg]. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
            The temperature [K].

        Returns
        -------
        float
            The hole concentration [cm<sup>&minus;3</sup>].
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nv(T) * fd1(-Ef / (k * T))

    def _intrinsic_charge_imbalance(self, Ef, T):
        """note: the function decreases monotonically with increasing Ef"""
        # This method is meant to be overridden.
        return self.p_concentration(Ef, T) - self.n_concentration(Ef, T)

    def _charge_imbalance(self, Ef, T):
        """note: the function decreases monotonically with increasing Ef"""
        return self.p_concentration(Ef, T) - self.n_concentration(Ef, T)

    def _solve_electroneutrality_equation(self, equation, T) -> float:
        eq = partial(equation, T=T)
        return bisect(eq, -self.Eg, 2 * self.Eg, xtol=1e-6 * self.Eg)  # noqa

    def intrinsic_fermi_level(self, T=300):
        """
        Determine the Fermi level from the condition of electroneutrality
        based on the intrinsic electron concentration.

        .. math::
            n_h - n_{e,intrinsic} = 0

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The Fermi level [erg].
        """
        return self._solve_electroneutrality_equation(self._intrinsic_charge_imbalance, T)

    def fermi_level(self, T=300):
        """
        Determine the Fermi level from the condition of electroneutrality:

        .. math::
            p - n = 0

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The Fermi level [erg].
        """
        return self._solve_electroneutrality_equation(self._charge_imbalance, T)

    def conductivity_type(self, *, T=None, Ef=None):
        """
        Tell the conductivity type.

        Returns
        -------
        str
            `'i'` -- the intrinsic type.
        """
        return 'i'


class DopedSemiconductor(Semiconductor):
    """
    A class to calculate properties of a doped semiconductor (extends `Semiconductor`).

    Attributes
    ----------
    Na : float
        The acceptor concentration [cm<sup>&minus;3</sup>].
    Ea : float
        The acceptor level [erg].
    Nd : float
        The donor concentration [cm<sup>&minus;3</sup>].
    Ed : float
        The donor level [erg].
    """

    def __init__(self, mat, Na, Ea, Nd, Ed):
        """
        Construct the necessary attributes for the `DopedSemiconductor` object.

        Parameters
        ----------
        mat : Semiconductor
            The base intrinsic (pure) semiconductor.
        Na : float
            The acceptor concentration [cm<sup>&minus;3</sup>].
        Ea : float
            The acceptor level [erg].
        Nd : float
            The donor concentration [cm<sup>&minus;3</sup>].
        Ed : float
            The donor level [erg].
        """
        super(DopedSemiconductor, self).__init__(mat.me, mat.mh, mat.Eg, mat.chi, mat.eps, mat.lattice)
        self.Na = Na
        self.Ea = Ea
        self.Nd = Nd
        self.Ed = Ed

    @staticmethod
    def from_materials(material, mobility, dopant, resistivity):
        # TODO: documentation
        N = concentration(resistivity, mobility)
        Na = Nd = 0
        Ea = 0
        Ed = material.Eg
        if dopant in material.acceptor_energy:
            Ea = material.acceptor_energy[dopant]
            Na = N
        elif dopant in material.donor_energy:
            Ed = material.donor_energy[dopant]
            Nd = N
        else:
            raise KeyError('Unknown dopant')
        return DopedSemiconductor(material, Na, Ea, Nd, Ed)

    def p_donor_concentration(self, Ef=None, T=300):
        r"""
        Calculate the concentration of positive donor ions.

        .. math::
            N_d^+ = N_d \cdot (1 - f(E_d))

        Parameters
        ----------
        Ef : float or None
            The Fermi level [erg]. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
            The temperature [K].

        Returns
        -------
        float
            The concentration of positive donor ions [cm<sup>&minus;3</sup>].
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Nd * (1 - fermi(self.Ed, Ef, T))

    def n_acceptor_concentration(self, Ef=None, T=300):
        r"""
        Calculate the concentration of negative acceptor ions.

        .. math::
            N_a^- = N_d f(E_a)

        Parameters
        ----------
        Ef : float or None
            The Fermi level [erg]. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
            The temperature [K].

        Returns
        -------
        float
            The concentration of negative acceptor ions [cm<sup>&minus;3</sup>].
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Na * fermi(self.Ea, Ef, T)

    def _charge_imbalance(self, Ef, T):
        """note: the function decreases monotonically with increasing Ef"""
        return self.p_concentration(Ef, T) + self.p_donor_concentration(Ef, T) \
               - self.n_concentration(Ef, T) - self.n_acceptor_concentration(Ef, T)

    def fermi_level(self, T=300):
        r"""
        Overrides `Semiconductor.fermi_level`

        Determine the Fermi level from the condition of electroneutrality:

        .. math::
            p + N_d^{+} - n - N_a^{-} = 0

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The Fermi level [erg].
        """
        return super().fermi_level(T)  # Note: the method _charge_imbalance() is overridden by DopedSemiconductor

    def conductivity_type(self, *, T=None, Ef=None):
        """
        Tell the conductivity type (overrides `Semiconductor.conductivity_type`).

        Parameters
        ----------
        Ef : float or None
            The Fermi level [erg]. If `None`, `Ef` is found via the `fermi_level` method.
        T : float or None
            The temperature [K]. If `None`, `T` is assigned 300 K.

        Returns
        -------
        str
            `'p'` -- the positive type; `'n'` -- the negative type.

        Raises
        ------
        ValueError
            Both `T` and `Ef` are specified (neither is `None`).
        """
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
            The work function [erg].
        """
        self.work_function = work_function


class ContactType(Enum):
    """
    An enumeration of contact types.
    """
    AUGMENTATION = 0,
    DEPLETION = 1,
    INVERSION = 2


class MSJunction:
    """
    A class to calculate properties of a contact between a metal and a semiconductor.

    Attributes
    ----------
    metal : Metal
        The metal.
    sc : Semiconductor
        The semiconductor.
    """

    def __init__(self, metal, sc):
        """
        Construct the necessary attributes for the `MSJunction` object.

        Parameters
        ----------
        metal : Metal
            The metal.
        sc : Semiconductor
            The semiconductor.
        """
        self.metal = metal
        self.sc = sc

    def delta_phi(self, T=300):
        r"""
        Calculate the difference between the exit potentials of the metal and the semiconductor (in units of voltage).

        .. math::
            \Delta \phi = - \frac{ E_g - E_f(T) + \chi - \Phi_M }{ e }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The difference of potentials [statV].
        """
        return -(self.sc.Eg - self.sc.fermi_level(T) + self.sc.chi - self.metal.work_function) / e

    def schottky_barrier(self):
        r"""
        Calculate the height of the Schottky barrier.

        .. math::
            \Phi_B = \frac{ \Phi_M - \chi }{ e }

        Returns
        -------
        float
            The height of the Schottky barrier [statV].
        """
        return (self.metal.work_function - self.sc.chi) / e

    def contact_type(self, T=300):
        """
        Determine the type of the contact.

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        ContactType
            The type of the contact.

        Raises
        ------
        NotImplementedError
            The semiconductor is of the intrinsic conductivity type.
        """
        ct = self.sc.conductivity_type(T=T)
        Ef = self.sc.fermi_level(T)
        dEf = -self.delta_phi(T) * e
        if ct == 'i':
            raise NotImplementedError('I do not know which cases are possible')
        else:
            if ct != self.sc.conductivity_type(Ef=Ef + dEf):
                return ContactType.INVERSION
            else:
                if ct == 'p':
                    return ContactType.AUGMENTATION if dEf < 0 else ContactType.DEPLETION
                if ct == 'n':
                    return ContactType.DEPLETION if dEf < 0 else ContactType.AUGMENTATION

    def full_depletion_width(self, T=300):
        r"""
        Calculate the width of the depletion region, using the approximation of full depletion.

        .. math::
            w = \sqrt{ \frac{ \epsilon \Delta \phi }{ 2 \pi e n } }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The width of the depletion region [cm].
        """
        return depletion_width(self.sc.eps, self.sc.n_concentration(T=T), self.delta_phi(T))

    def debye_length(self, T=300):
        r"""
        Calculate the Debye length (screening length).

        .. math::
            \lambda_D = \sqrt{ \frac{ \epsilon k T }{ 4 \pi e^2 n_e } }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The Debye length [cm].
        """
        return debye_length(self.sc.eps, self.sc.n_concentration(T=T), T)


class PNJunction:
    """
    A class to calculate properties of a p-n junction.

    Attributes
    ----------
    mat : Semiconductor
        The base intrinsic (pure) semiconductor.
    n_mat : DopedSemiconductor
        The n-type doped semiconductor.
    p_mat : DopedSemiconductor
        The p-type doped semiconductor.
    """

    def __init__(self, mat, Na, Ea, Nd, Ed):
        """
        Construct the necessary attributes for the `PNJunction` object.

        Parameters
        ----------
        mat : Semiconductor
            The base intrinsic (pure) semiconductor.
        Na : float
            The acceptor concentration [cm<sup>&minus;3</sup>].
        Ea : float
            The acceptor level [erg].
        Nd : float
            The donor concentration [cm<sup>&minus;3</sup>].
        Ed : float
            The donor level [erg].
        """
        if Ea is None:
            Ea = 0
        if Ed is None:
            Ed = mat.Eg
        self.mat = mat
        self.p_mat = DopedSemiconductor(mat, Na, Ea, 0, Ed)
        self.n_mat = DopedSemiconductor(mat, 0, Ea, Nd, Ed)

    def delta_phi(self, T=300):
        r"""
        Calculate the difference between the donor and acceptor Fermi potentials.

        .. math::
            \Delta \phi = \frac{ E_{f,n}(T) - E_{f,p}(T) }{ e }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The difference of potentials [statV].
        """
        return (self.n_mat.fermi_level(T) - self.p_mat.fermi_level(T)) / e


class PNJunctionNonDegenerate(PNJunction):
    """
    A class to calculate properties of a p-n junction in the non-degenerate semiconductor approximation
    (extends `PNJunction`).
    """

    def pn(self, voltage, T=300):
        r"""
        Calculate the product of the hole and electron concentrations.

        .. math::
            p \cdot n = n_i^2 \exp{\frac{e V}{k T}}

        Parameters
        ----------
        voltage : float
            The voltage at the p-n junction [statV].
        T : float
            The temperature [K].

        Returns
        -------
        float
            The hole concentration multiplied by the electron concentration [cm<sup>&minus;6</sup>].
        """
        return self.mat.i_concentration(T) ** 2 * exp(e * voltage / (k * T))

    def n_p(self, voltage, T=300):
        r"""
        Calculate the electron concentration in the p-semiconductor.

        .. math::
            n_p = \frac{n_i^2}{p_p} \exp{\frac{e V}{k T}}

        Parameters
        ----------
        voltage : float
            The voltage at the p-n junction [statV].
        T : float
            The temperature [K].

        Returns
        -------
        float
            The electron concentration in the p-semiconductor [cm<sup>&minus;3</sup>].
        """
        return self.pn(voltage, T) / self.p_mat.p_concentration(T=T)

    def p_n(self, voltage, T=300):
        r"""
        Calculate the hole concentration in the n-semiconductor.

        .. math::
            p_n = \frac{n_i^2}{n_n} \exp{\frac{e V}{k T}}

        Parameters
        ----------
        voltage : float
            The voltage at the p-n junction [statV].
        T : float
            The temperature [K].

        Returns
        -------
        float
            The hole concentration in the n-semiconductor [cm<sup>&minus;3</sup>].
        """
        return self.pn(voltage, T) / self.n_mat.n_concentration(T=T)

    def j0_p(self, diffusivity, diffusion_length):
        r"""
        Calculate the hole current density without external voltage (the dark current).

        .. math::
            J_{0p} = \frac{e D_p p_{n0}}{L_p}

        Parameters
        ----------
        diffusivity : float
            The hole diffusivity [cm<sup>2</sup> s<sup>&minus;1</sup>].
        diffusion_length : float
            The hole diffusion length [cm].

        Returns
        -------
        float
            The hole current density [statA cm<sup>&minus;2</sup>].
        """
        return e * diffusivity * self.p_n(0) / diffusion_length

    def j0_n(self, diffusivity, diffusion_length):
        r"""
        Calculate the electron current density without external voltage (the dark current).

        .. math::
            J_{0n} = \frac{e D_n n_{p0}}{L_n}

        Parameters
        ----------
        diffusivity : float
            The electron diffusivity [cm<sup>2</sup> s<sup>&minus;1</sup>].
        diffusion_length : float
            The electron diffusion length [cm].

        Returns
        -------
        float
            The electron current density [statA cm<sup>&minus;2</sup>].
        """
        return e * diffusivity * self.p_n(0) / diffusion_length

    def current_p(self, diffusivity, diffusion_length, voltage, T=300):
        r"""
        Calculate the hole current density.

        .. math::
            J_p = J_{0p}\left[\exp{\frac{e V}{k T}} - 1\right]
                = \frac{e D_p p_{n0}}{L_p}\left[\exp{\frac{e V}{k T}} - 1\right]

        Parameters
        ----------
        diffusivity : float
            The hole diffusivity [cm<sup>2</sup> s<sup>&minus;1</sup>].
        diffusion_length : float
            The hole diffusion length [cm].
        voltage : float
            The voltage at the p-n junction [statV].
        T : float
            The temperature [K].

        Returns
        -------
        float
            The hole current density [statA cm<sup>&minus;2</sup>].
        """
        return self.j0_p(diffusivity, diffusion_length) * (exp(e * voltage / (k * T)) - 1)

    def current_n(self, diffusivity, diffusion_length, voltage, T=300):
        r"""
        Calculate the electron current density.

        .. math::
            J_n = J_{0n}\left[\exp{\frac{e V}{k T}} - 1\right]
                = \frac{e D_n n_{p0}}{L_n}\left[\exp{\frac{e V}{k T}} - 1\right]

        Parameters
        ----------
        diffusivity : float
            The electron diffusivity [cm<sup>2</sup> s<sup>&minus;1</sup>].
        diffusion_length : float
            The electron diffusion length [cm].
        voltage : float
            The voltage at the p-n junction [statV].
        T : float
            The temperature [K].

        Returns
        -------
        float
            The electron current density [statA cm<sup>&minus;2</sup>].
        """
        return self.j0_n(diffusivity, diffusion_length) * (exp(e * voltage / (k * T)) - 1)


class PNJunctionFullDepletion(PNJunction):
    """
    A class to calculate properties of a p-n junction exhibiting full depletion (extends `PNJunction`).
    """

    def __init__(self, mat, Na, Ea, Nd, Ed):
        """
        Construct the necessary attributes for the `PNJunctionFullDepletion` object
        (calls the `PNJunction` constructor).
        """
        super().__init__(mat, Na, Ea, Nd, Ed)

    def _df_Na_Nd(self, T):
        return self.delta_phi(T), \
               self.p_mat.n_acceptor_concentration(T=T), \
               self.n_mat.p_donor_concentration(T=T)

    def _w_tmp(self, T):
        df, a, d = self._df_Na_Nd(T)
        return self.p_mat.eps / (2 * pi * e) * df, a, d

    def delta_phi_n(self, T=300):
        r"""
        Calculate the difference of potentials in the n-type semiconductor.

        .. math::
            \Delta\phi_n = \Delta\phi \frac{ N_a^- }{ N_a^- + N_d^+ }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The difference of potentials [statV].
        """
        df, a, d = self._df_Na_Nd(T)
        return df * a / (a + d)

    def delta_phi_p(self, T=300):
        r"""
        Calculate the difference of potentials in the p-type semiconductor.

        .. math::
            \Delta\phi_p = \Delta\phi \frac{ N_d^+ }{ N_a^- + N_d^+ }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The difference of potentials [statV].
        """
        df, a, d = self._df_Na_Nd(T)
        return df * d / (a + d)

    def w(self, T=300):
        r"""
        Calculate the full depletion width.

        .. math::
            w = \sqrt{ \frac{ \epsilon }{ 2 \pi e } \Delta\phi \frac{ N_a^- + N_d^+ }{ N_a^-  N_d^+ } }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The full depletion width [cm].
        """
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * (a + d) / (a * d))

    def w_n(self, T=300):
        r"""
        Calculate the width of the depletion layer inside the n-type semiconductor.

        .. math::
            w = \sqrt{ \frac{ \epsilon }{ 2 \pi e } \Delta\phi \frac{ N_a^- }{ N_d^+ } \frac{ 1 }{ N_a^-  N_d^+ } }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The width of the depletion layer inside the n-type semiconductor [cm].
        """
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * a / d / (a + d))

    def w_p(self, T=300):
        r"""
        Calculate the width of the depletion layer inside the p-type semiconductor.

        .. math::
            w = \sqrt{ \frac{ \epsilon }{ 2 \pi e } \Delta\phi \frac{ N_d^+ }{ N_a^- } \frac{ 1 }{ N_a^-  N_d^+ } }

        Parameters
        ----------
        T : float
            The temperature [K].

        Returns
        -------
        float
            The width of the depletion layer inside the p-type semiconductor [cm].
        """
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * d / a / (a + d))


class PeriodicPotentialModel(ABC):
    """"""

    # TODO: add documentation
    def __init__(self, u_min, period):
        self.u_min = u_min
        self.period = period

    def equation(self, energy, k, m):
        return self.equation_left_part(energy, m) - self.equation_right_part(k)

    def equation_right_part(self, k):
        return cos(k * self.period)

    @abstractmethod
    def equation_left_part(self, energy, m):
        """The part of the equation independent of k"""

    def find_lower_band_range(self, m, xtol_coarse, xtol_fine):
        assert xtol_coarse > 0 and xtol_fine > 0

        # Find band with minimal energy (starts search from self.u_min)
        def f1(energy):
            return self.equation_left_part(energy, m) - 1

        def f2(energy):
            return self.equation_left_part(energy, m) + 1

        e_start = find_nth_function_zero(f1, self.u_min, xtol_coarse, xtol_fine, num=0)
        e_end = find_nth_function_zero(f2, e_start - xtol_coarse / 2, xtol_coarse, xtol_fine, num=0)

        return e_start, e_end

    def get_energy(self, k, m, bracket, xtol):
        start, stop = bracket
        return bisect(self.equation, start, stop, args=(k, m), xtol=xtol)  # noqa

    def get_k(self, energy, m):
        # Returns k multiplied by period
        val = self.equation_left_part(energy, m)
        if not -1 <= val <= 1:
            return None
        return math.acos(val)

    def get_ks(self, es, m):  # vectorized version of get_k
        # Returns k multiplied by period
        vs = np.vectorize(self.equation_left_part)(es, m)
        return np.arccos(np.clip(vs, -1, 1))


class KronigPenneyModel(PeriodicPotentialModel):
    """"""

    # TODO: add documentation (add equations to the description of the class)
    def __init__(self, a, b, u0):
        assert a > 0 and b > 0
        super().__init__(min(0, u0), a + b)
        self.a = a
        self.b = b
        self.u0 = u0

    def equation(self, energy, k, m):
        r"""
        a - is the width of area where potential energy is U0

        b - is the width of area where potential energy is 0

        .. math::
             cos(\alpha a) cos(\beta b) -
            \frac{\alpha^2 + \beta^2}{2 \alpha \beta} sin(\alpha a) sin(\beta b) = cos(k (a + b))
        .. math::
            \alpha^2 = \frac{2 m (E - U_0)}{\hbar^2}; \beta^2 = \frac{2 m E}{\hbar^2}
        """
        return super(KronigPenneyModel, self).equation(energy, k, m)

    def equation_left_part(self, energy, m):
        alf = cmath.sqrt(2 * m / h_bar ** 2 * (energy - self.u0))
        bet = cmath.sqrt(2 * m / h_bar ** 2 * energy)
        first = cmath.cos(alf * self.a) * cmath.cos(bet * self.b)
        if bet == 0:
            second = - alf / 2 * self.b * cmath.sin(alf * self.a)
        elif alf == 0:
            second = - bet / 2 * self.a * cmath.sin(bet * self.b)
        else:
            second = - (alf ** 2 + bet ** 2) / (2 * alf * bet) * cmath.sin(alf * self.a) * cmath.sin(bet * self.b)
        return (first + second).real


class DiracCombModel(PeriodicPotentialModel):
    """"""

    # TODO: add documentation
    def __init__(self, a, G):
        assert a > 0
        super().__init__(-1e-10 * eV, a)
        self.G = G

    def equation(self, energy, k, m):
        r"""
        .. math::
            cos(\beta b) + \sqrt{\frac{m}{2 \hbar^2 E}} G sin(\beta b) = cos(k b)
        .. math::
            \beta^2 = \frac{2 m E}{\hbar^2}

        Alternative formula used in our course previously (note: different notation for alpha and beta)

        .. math::
            cos(\alpha a) + \frac{2 m G}{k \hbar^2} sin(\alpha a) = cos(k a)
        .. math::
            \alpha^2 = \frac{2 m E}{\hbar^2}
        """
        return super(DiracCombModel, self).equation(energy, k, m)

    def equation_left_part(self, energy, m):
        bet = cmath.sqrt(2 * m / h_bar ** 2 * energy)
        first = cmath.cos(bet * self.period)
        if energy == 0:
            second = m * self.G * self.period / h_bar ** 2
        else:
            second = cmath.sqrt(m / (2 * h_bar ** 2 * energy)) * self.G * cmath.sin(bet * self.period)
        return (first + second).real


class JFET:
    def __init__(self, material, Nd, mobility, a, b, L):
        self.material = material
        self.Nd = Nd
        self.mobility = mobility
        self.a = a
        self.b = b
        self.L = L

    def Vp(self):
        r"""
        .. math::
            V_p = \frac{2 \pi e N_d a^3}{\epsilon}
        """
        return 2 * pi * e * self.Nd * self.a ** 2 / self.material.eps

    def Ip(self):
        r"""
        .. math::
            I_p = \frac{4 \pi e^2 \mu N_d^2 a^3 b}{3 \epsilon L}
        """
        return 4 * pi * e ** 2 * self.mobility * self.b * self.Nd ** 2 * self.a ** 3 / (3 * self.material.eps * self.L)

    def Id_sat(self, Vg, delta_phi):
        r"""
        .. math::
            I_{D,sat} = I_P ( 1 - 3 \frac{Vg + \Delta \phi}{V_p} + 2 \frac{ ( V_g + \Delta \phi )^{\frac{3}{2}}}{V_p^{\frac{3}{2}}})
        """
        Vp = self.Vp()
        Ip = self.Ip()
        v_sum = Vg + delta_phi
        return Ip * (1 - 3 * v_sum / Vp + 2 * (v_sum / Vp) ** (3 / 2))
