"""
This module contains useful functions.
"""

from math import exp

import numpy as np

from fompy.constants import k
from fompy.util import fermi_dirac

fd1 = np.vectorize(fermi_dirac.fd1)


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
