from enum import Enum
from functools import partial
from math import pi, sqrt

from scipy.optimize import bisect

from fompy.constants import e, k, h_bar
from fompy.functions import fermi, fd1


def conductivity(n, n_mob, p, p_mob):
    return e * (n * n_mob + p * p_mob)


# w depletion area width
def depletion_width(eps, n, d_phi):
    """Ширина зоны между металлом и полупроводником"""
    return sqrt(eps * d_phi / (2 * pi * e * n))


# L_D screening length
def debye_length(eps, n, T):
    """Длина Дебая"""
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
        """
        Used to calculate Nc and Nv

        m_eff -- Effective mass of density of states
        T -- Temperature
        """
        return 2 * (2 * pi * m_eff * k * T / (2 * pi * h_bar) ** 2) ** (3 / 2)

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
        return Semiconductor.effective_state_density(self.me, T)

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
        return Semiconductor.effective_state_density(self.mh, T)

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
        super(DopedSemiconductor, self).__init__(mat.me, mat.mh, mat.Eg, mat.chi, mat.eps)
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
        return self.Nd * (1 - fermi(self.Ed, Ef, T))

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
        return self.Na * fermi(self.Ea, Ef, T)

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


class ContactType(Enum):
    AUGMENTATION = 0,
    DEPLETION = 1,
    INVERSION = 2


class MetalSemiconductorContact:
    def __init__(self, metal, sc):
        self.metal = metal
        self.sc = sc

    def delta_phi(self, T=300):
        """Returns difference of exit potential of metal and semiconductor (in units of voltage)"""
        return -(self.sc.Eg - self.sc.fermi_level(T) + self.sc.chi - self.metal.work_function) / e

    def contact_type(self, T=300):
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

    def __init__(self, mat, Na, Ea, Nd, Ed):
        if Ea is None:
            Ea = 0
        if Ed is None:
            Ed = mat.Eg
        self.mat = mat
        self.p_mat = DopedSemiconductor(mat, Na, Ea, 0, Ed)
        self.n_mat = DopedSemiconductor(mat, 0, Ea, Nd, Ed)

    def delta_phi(self, T=300):
        return (self.n_mat.fermi_level(T) - self.p_mat.fermi_level(T)) / e


class PNJunctionFullDepletion(PNJunction):

    def __init__(self, mat, Na, Ea, Nd, Ed):
        super().__init__(mat, Na, Ea, Nd, Ed)

    def _df_Na_Nd(self, T):
        return self.delta_phi(T), \
               self.p_mat.n_acceptor_concentration(T=T), \
               self.n_mat.p_donor_concentration(T=T)

    def _w_tmp(self, T):
        df, a, d = self._df_Na_Nd(T)
        return self.p_mat.eps / (2 * pi * e) * df, a, d

    def delta_phi_n(self, T=300):
        df, a, d = self._df_Na_Nd(T)
        return df * a / (a + d)

    def delta_phi_p(self, T=300):
        df, a, d = self._df_Na_Nd(T)
        return df * d / (a + d)

    def w(self, T=300):
        r"""
        Compute the full depletion width.

        .. math::
            w = \sqrt{ { \epsilon \over 2 \pi e} \Delta\phi { N_a^- + N_d^+  \over  N_a^-  N_d^+ } }

        Parameters
        ----------
        T=300 : float
            The temperature.

        Returns
        -------
        float
            The full depletion width.

        """
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * (a + d) / (a * d))

    def w_n(self, T=300):
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * a / d / (a + d))

    def w_p(self, T=300):
        tmp, a, d = self._w_tmp(T)
        return sqrt(tmp * d / a / (a + d))
