"""
This module contains instances of classes from `fompy.models` for particular materials.
"""

from fompy.constants import me, eV, angstrom, amu
from fompy.models import Semiconductor, DiamondLikeLattice

# http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html
Si = Semiconductor(
    0.36 * me, 0.81 * me, 1.12 * eV, 4.05 * eV, 11.7,
    DiamondLikeLattice(5.4307 * angstrom, 28 * amu))
"""A `Semiconductor` object for silicon at 300 K."""

# http://www.ioffe.ru/SVA/NSM/Semicond/Ge/index.html
Ge = Semiconductor(
    0.22 * me, 0.34 * me, 0.661 * eV, 4.0 * eV, 16.2,
    DiamondLikeLattice(5.660 * angstrom, 72.6 * amu))
"""A `Semiconductor` object for germanium at 300 K."""

# http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/index.html
GaAs = Semiconductor(
    0.063 * me, 0.53 * me, 1.424 * eV, 4.07 * eV, 12.9,
    DiamondLikeLattice(5.6533 * angstrom, (69.723 + 74.922) / 2 * amu))  # Gamma-valley
"""A `Semiconductor` object for gallium arsenide at 300 K in Gamma valley."""

# TODO: add more materials
