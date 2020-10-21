"""
This module contains useful functions.
"""

from math import exp

from fompy.constants import k
from fompy.util.fermi_dirac import fd1  # noqa


def fermi(E, Ef, T):
    r"""
    Calculate the Fermi-Dirac distribution.

    .. math::
        f(E) = \frac{ 1 }{ 1 + \exp\left( \frac{ E - E_f }{ k T } \right) }

    Parameters
    ----------
    E : float
        The energy level.
    Ef : float
        The Fermi level.
    T : float
        The temperature.

    Returns
    -------
    float
        The Fermi-Dirac distribution of `E`.
    """
    return 1 / (1 + exp((E - Ef) / (k * T)))
