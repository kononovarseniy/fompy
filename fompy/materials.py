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

from fompy.constants import me, eV
from fompy.models import Semiconductor

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
