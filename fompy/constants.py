"""
This file contains values of frequently used physical constants.

All values are given in Gaussian unit system (cgs units), rounded to four
significant digits according to the rules of arithmetic.
"""

################################################################################

c = 2.998e10
"""The Speed of light in vacuum"""

################################################################################

k = 1.381e-16
"""The Boltzmann constant"""

Na = 6.022e23
"""The Avogadro constant"""

R = 8.314e7
"""The [universal] gas constant"""

sigma = 5.670e-5
"""The Stefan-Boltzmann constant"""

################################################################################

me = 9.109e-28
"""Electron mass"""

mp = 1.673e-24
"""Proton mass"""

mn = 1.675e-24
"""Neutron mass"""

e = 4.803e-10
"""Electron charge"""

h = 6.626e-27
"""The Planck constant"""

h_bar = 1.055e-27
"""Reduced Planck constant"""

################################################################################

eV = 1.602e-12
"""One electronvolt expressed in ergs"""

eV_m = 1.783e-33
"""Off-system unit of mass equal to eV / c^2"""

eV_T = 1.160e4
"""Off-system unit of temperature equal to eV / k"""

amu = 1.661e-24
"""Atomic Mass Unit (Dalton)"""

angstrom = 1e-8
"""One angstrom expressed in centimeters"""

################################################################################

volt = 1e8 / c
"""One volt"""

ampere = 1e-1 * c
"""One ampere"""

ohm = 1e9 / c**2
"""One ohm"""

farad = 1e-9 * c**2
"""One farad"""

henry = 1e9 / c**2
"""One henry"""

################################################################################

Ry = e**4*me/(2*h_bar**2)

a0 = h_bar**2/(me*e**2)

