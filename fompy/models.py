"""
This module contains classes useful for calculating properties of different materials.

Notes
-----
All energies are counted from the valence band edge Ev.
"""

from enum import Enum
from functools import partial
from math import pi, sqrt

from scipy.optimize import bisect

from fompy.constants import e, k, h_bar
from fompy.functions import fermi, fd1


def conductivity(n, n_mob, p, p_mob):
    r"""
    Calculate the conductivity of a material.

    .. math::
        \sigma = e (n_n \mu_n + n_p \mu_p)

    Parameters
    ----------
    n : float
        The concentration of electrons.
    n_mob : float
        The electron mobility.
    p : float
        The concentration of holes.
    p_mob : float
        The hole mobility.

    Returns
    -------
    float
        The conductivity.
    """
    return e * (n * n_mob + p * p_mob)


# w depletion area width
def depletion_width(eps, n, d_phi):
    r"""
    Calculate the width of the depletion region.

    .. math::
        w = \sqrt{ \frac{ \epsilon \Delta \phi }{ 2 \pi e n } }

    Parameters
    ----------
    eps : float
        The dielectric constant.
    n : float
        The concentration of charge carriers.
    d_phi : float
        The difference of potentials.

    Returns
    -------
    float
        The width of the depletion region.
    """
    return sqrt(eps * d_phi / (2 * pi * e * n))


# L_D screening length
def debye_length(eps, n, T):
    r"""
    Calculate the Debye length.

    .. math::
        \lambda_D = \sqrt{ \frac{ \epsilon k T }{ 4 \pi e^2 n } }

    Parameters
    ----------
    eps : float
        The dielectric constant.
    n : float
        The concentration of charge carriers.
    T : float
        The temperature.

    Returns
    -------
    float
        The Debye length.
    """
    return sqrt(eps * k * T / (4 * pi * e ** 2 * n))


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
    eps : float
        The dielectric constant.
    """

    def __init__(self, me_eff, mh_eff, Eg, chi, eps):
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
        eps : float
            The dielectric constant.
        """
        self.me = me_eff
        self.mh = mh_eff
        self.Eg = Eg
        self.chi = chi
        self.eps = eps

    @staticmethod
    def effective_state_density(m_eff, T):
        r"""
        Calculate the effective density of states.

        .. math::
            N = 2 \left( \frac{2 \pi m_{eff} k T }{ (2 \pi \hbar)^2 } \right)^{3/2}

        Parameters
        ----------
        m_eff : float
            The effective mass.
        T : float
            The temperature.

        Returns
        -------
        float
            The effective density of states.
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
            The temperature.

        Returns
        -------
        float
            The effective density of states for electrons in the conduction band.
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
            The temperature.

        Returns
        -------
        float
            The effective density of states for holes in the valence band.
        """
        return Semiconductor.effective_state_density(self.mh, T)

    def n_concentration(self, Ef=None, T=300):
        r"""
        Calculate the electron concentration.

        .. math::
            n_n = N_c(T) \Phi_{1/2}\left( \frac{ E_f - E_g }{ k T } \right)

        Parameters
        ----------
        Ef : float or None
            The Fermi level. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
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
        r"""
        Calculate the hole concentration.

        .. math::
            n_p = N_v(T) \Phi_{1/2}\left( \frac{ - E_f }{ k T } \right)

        Parameters
        ----------
        Ef : float or None
            The Fermi level. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
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

        .. math::
            n_p - n_n = 0

        Parameters
        ----------
        T : float
            The temperature.

        Returns
        -------
        float
            The Fermi level.
        """
        # TODO: Not sure if it works correctly with Nd!=0 and Na!=0
        return bisect(partial(self._charge_imbalance, T=T), 0, self.Eg, xtol=1e-6 * self.Eg)

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
            The base intrinsic (pure) semiconductor.
        Na : float
            The acceptor concentration.
        Ea : float
            The acceptor level.
        Nd : float
            The donor concentration.
        Ed : float
            The donor level.
        """
        super(DopedSemiconductor, self).__init__(mat.me, mat.mh, mat.Eg, mat.chi, mat.eps)
        self.Na = Na
        self.Ea = Ea
        self.Nd = Nd
        self.Ed = Ed

    def p_donor_concentration(self, Ef=None, T=300):
        r"""
        Calculate the concentration of positive donor ions.

        .. math::
            N_d^+ = N_d \cdot (1 - f(E_d))

        Parameters
        ----------
        Ef : float or None
            The Fermi level. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
            The temperature.

        Returns
        -------
        float
            The concentration of positive donor ions.
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
            The Fermi level. If `None`, `Ef` is found via the `fermi_level` method.
        T : float
            The temperature.

        Returns
        -------
        float
            The concentration of negative acceptor ions.
        """
        if Ef is None:
            Ef = self.fermi_level(T)
        return self.Na * fermi(self.Ea, Ef, T)

    def _charge_imbalance(self, Ef, T):
        return self.p_concentration(Ef, T) + self.p_donor_concentration(Ef, T) \
               - self.n_concentration(Ef, T) - self.n_acceptor_concentration(Ef, T)

    def conductivity_type(self, *, T=None, Ef=None):
        """
        Tell the conductivity type (overrides `Semiconductor.conductivity_type`).

        Parameters
        ----------
        Ef : float or None
            The Fermi level. If `None`, `Ef` is found via the `fermi_level` method.
        T : float or None
            The temperature. If `None`, `T` is assigned 300 K.

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
            The work function.
        """
        self.work_function = work_function


class ContactType(Enum):
    """
    An enumeration of contact types.
    """
    AUGMENTATION = 0,
    DEPLETION = 1,
    INVERSION = 2


class MetalSemiconductorContact:
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
        Construct the necessary attributes for the `MetalSemiconductorContact` object.

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
            \Delta \phi = - \frac{ E_g - E_f(T) + \chi - W }{ e }

        Parameters
        ----------
        T : float
            The temperature.

        Returns
        -------
        float
            The difference of potentials.
        """
        return -(self.sc.Eg - self.sc.fermi_level(T) + self.sc.chi - self.metal.work_function) / e

    def contact_type(self, T=300):
        """
        Determine the type of the contact.

        Parameters
        ----------
        T : float
            The temperature.

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
            The acceptor concentration.
        Ea : float
            The acceptor level.
        Nd : float
            The donor concentration.
        Ed : float
            The donor level.
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
            The temperature.

        Returns
        -------
        float
            The difference of potentials.
        """
        return (self.n_mat.fermi_level(T) - self.p_mat.fermi_level(T)) / e


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
        Calculate the difference of potentials for the negative depletion layer.

        .. math::
            \Delta\phi_n = \Delta\phi \frac{ N_a^- }{ N_a^- + N_d^+ }

        Parameters
        ----------
        T : float
            The temperature.

        Returns
        -------
        float
            The difference of potentials.
        """
        df, a, d = self._df_Na_Nd(T)
        return df * a / (a + d)

    def delta_phi_p(self, T=300):
        r"""
        Calculate the difference of potentials for the positive depletion layer.

        .. math::
            \Delta\phi_p = \Delta\phi \frac{ N_d^+ }{ N_a^- + N_d^+ }

        Parameters
        ----------
        T : float
            The temperature.

        Returns
        -------
        float
            The difference of potentials.
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
            The temperature.

        Returns
        -------
        float
            The full depletion width.
        """
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * (a + d) / (a * d))

    def w_n(self, T=300):
        r"""
        Calculate the width of the negative depletion layer.

        .. math::
            w = \sqrt{ \frac{ \epsilon }{ 2 \pi e } \Delta\phi \frac{ N_a^- }{ N_d^+ } \frac{ 1 }{ N_a^-  N_d^+ } }

        Parameters
        ----------
        T : float
            The temperature.

        Returns
        -------
        float
            The width of the negative depletion layer.
        """
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * a / d / (a + d))

    def w_p(self, T=300):
        r"""
        Calculate the width of the negative depletion layer.

        .. math::
            w = \sqrt{ \frac{ \epsilon }{ 2 \pi e } \Delta\phi \frac{ N_d^+ }{ N_a^- } \frac{ 1 }{ N_a^-  N_d^+ } }

        Parameters
        ----------
        T : float
            The temperature.

        Returns
        -------
        float
            The width of the negative depletion layer.
        """
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * d / a / (a + d))
