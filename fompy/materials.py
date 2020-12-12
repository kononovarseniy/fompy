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
Si.donor_energy['As'] = Si.Eg - 0.054 * eV
Si.donor_energy['P'] = Si.Eg - 0.045 * eV
Si.donor_energy['Sb'] = Si.Eg - 0.043 * eV
Si.acceptor_energy['Al'] = 0.072 * eV
Si.acceptor_energy['B'] = 0.045 * eV
Si.acceptor_energy['Ga'] = 0.074 * eV
Si.acceptor_energy['In'] = 0.157 * eV

# http://www.ioffe.ru/SVA/NSM/Semicond/Ge/index.html
Ge = Semiconductor(
    0.22 * me, 0.34 * me, 0.661 * eV, 4.0 * eV, 16.2,
    DiamondLikeLattice(5.660 * angstrom, 72.6 * amu))
"""A `Semiconductor` object for germanium at 300 K."""
Ge.donor_energy['As'] = Ge.Eg - 0.014 * eV
Ge.donor_energy['P'] = Ge.Eg - 0.013 * eV
Ge.donor_energy['Sb'] = Ge.Eg - 0.010 * eV
Ge.donor_energy['Bi'] = Ge.Eg - 0.013 * eV
Ge.donor_energy['Li'] = Ge.Eg - 0.093 * eV
Ge.acceptor_energy['Al'] = 0.011 * eV
Ge.acceptor_energy['B'] = 0.011 * eV
Ge.acceptor_energy['Ga'] = 0.011 * eV
Ge.acceptor_energy['In'] = 0.012 * eV
Ge.acceptor_energy['Tl'] = 0.013 * eV

# http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/index.html
GaAs = Semiconductor(
    0.063 * me, 0.53 * me, 1.424 * eV, 4.07 * eV, 12.9,
    DiamondLikeLattice(5.6533 * angstrom, (69.723 + 74.922) / 2 * amu))  # Gamma-valley
"""A `Semiconductor` object for gallium arsenide at 300 K in Gamma valley."""
GaAs.donor_energy['S'] = GaAs.Eg - 0.006 * eV
GaAs.donor_energy['Se'] = GaAs.Eg - 0.006 * eV
GaAs.donor_energy['Si'] = GaAs.Eg - 0.006 * eV
GaAs.donor_energy['Ge'] = GaAs.Eg - 0.006 * eV
GaAs.donor_energy['Sn'] = GaAs.Eg - 0.006 * eV
GaAs.donor_energy['Te'] = GaAs.Eg - 0.03 * eV
GaAs.donor_energy['C'] = 0.02 * eV
GaAs.donor_energy['Ge'] = 0.03 * eV
GaAs.donor_energy['Zn'] = 0.025 * eV
GaAs.donor_energy['Sn'] = 0.02 * eV

# TODO: add more materials
