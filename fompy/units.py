"""
This module contains classes for parsing units of measurement.
"""

from fractions import Fraction

from fompy import constants
from fompy.util.parser import Tokenizer, TokenList

TOKENIZER = Tokenizer({
    'unit': r'[a-zA-Z]+',
    'num': r'([0-9]+)|([\-+][0-9]+)',
    'div': r'/',
    'optional': r'([*^])|(\s+)',
})


class SimpleUnit:
    """
    A class to convert a simple unit of measurement into cgs.

    Attributes
    ----------
    name : str
    prefix : str or None
    multiplier : float
        The multiplier associated with the prefix.
    value : float
        The value of the unit in cgs.
    power : int
    """

    def __init__(self, name, power):
        """
        Construct the necessary attributes for the `SimpleUnit` object.

        Parameters
        ----------
        name : str
        power : int
        """
        prefix = None
        multiplier = 1
        unit_name = name
        for p, n in PREFIXES.items():
            if len(p) < len(name) and name.startswith(p):
                if prefix is None or len(p) > len(prefix):
                    unit_name = name[len(p):]
                    prefix = p
                    multiplier = n

        self.name = unit_name
        self.prefix = prefix
        self.multiplier = multiplier
        self.value = UNITS[unit_name]
        self.power = power

    def get_number(self):
        """
        Get the value of the unit expressed in cgs.

        Returns
        -------
        float
        """
        return (self.multiplier * self.value) ** self.power

    def __str__(self):
        res = (self.prefix or '') + self.name
        if self.power != 1:
            res += '^' + str(self.power)
        return res


class CompositeUnit:
    """
    A class to convert a composite unit of measurement into cgs.

    Attributes
    ----------
    nom : iterable(CompositeUnit, SimpleUnit)
        The numerator.
    denom : iterable(CompositeUnit, SimpleUnit)
        The denominator.
    """

    def __init__(self, nom, denom):
        self.nom = nom
        self.denom = denom

    def get_number(self):
        """
        Get the value of the unit expressed in cgs.

        Returns
        -------
        float
        """
        res = 1
        for u in self.nom:
            res *= u.get_number()
        for u in self.denom:
            res /= u.get_number()
        return res

    def __str__(self):
        if len(self.nom) > 0:
            res = ' '.join(map(str, self.nom))
        else:
            res = '1'
        if len(self.denom) > 0:
            res += ' / ' + ' '.join(map(str, self.denom))
        return res


def parse_unit(text):
    """
    Extract a unit of measurement out of `text`.

    Parameters
    ----------
    text : str

    Returns
    -------
    CompositeUnit

    Raises
    ------
    SyntaxError
        The parser has encountered an unexpected token.
    """
    ts = TokenList(filter(lambda t: t[0] != 'optional', TOKENIZER.tokenize(text)))
    nom_list = []
    denom_list = []
    cur_list = nom_list
    while not ts.eof:
        tn, tv = ts.get()
        if tn == 'unit':
            unit_name = tv
            pow_nom = 1
            pow_denom = 1

            tn, tv = ts.get(1)
            if tn == 'num':
                pow_nom = int(tv)
                t1n, t1v = ts.get(2)
                t2n, t2v = ts.get(3)
                if t1n == 'div' and t2n == 'num':
                    pow_denom = int(t2v)
                    ts.position += 4
                else:
                    ts.position += 2
            else:
                ts.position += 1
            cur_list.append(SimpleUnit(unit_name, Fraction(pow_nom, pow_denom)))
        elif tn == 'num' and int(tv) == 1:  # Allow writing things like 1 / s
            ts.position += 1
        elif tn == 'div' and cur_list == nom_list:
            cur_list = denom_list
            ts.position += 1
        else:
            raise SyntaxError(f'Unexpected token {tn} {tv}')
    return CompositeUnit(nom_list, denom_list)


def unit(text):
    """
    Get the value of a unit of measurement in `text`.

    Parameters
    ----------
    text : str

    Returns
    -------
    float
        The value of the unit of measurement expressed in cgs.
    """
    return parse_unit(text).get_number()


def to_unit(val, unit_text):
    # TODO: documentation
    return val / unit(unit_text)


def from_unit(val, unit_text):
    # TODO: documentation
    return val * unit(unit_text)


def _register_unit(name, text):
    UNITS[name] = unit(text)


PREFIXES = {
    'y': 1e-24,
    'z': 1e-21,
    'a': 1e-18,
    'f': 1e-15,
    'p': 1e-12,
    'n': 1e-9,
    'u': 1e-6,
    'm': 1e-3,
    'c': 1e-2,
    'd': 1e-1,
    'da': 1e1,
    'h': 1e2,
    'k': 1e3,
    'M': 1e6,
    'G': 1e9,
    'T': 1e12,
    'P': 1e15,
    'E': 1e18,
    'Z': 1e21,
    'Y': 1e24
}

UNITS = {
    'm': 1e2,
    's': 1,
    'g': 1,
    'K': 1,
    'A': constants.ampere,

    'eV': constants.eV,
    'eV_m': constants.eV_m,
    'eV_T': constants.eV_T,
}

_register_unit('Hz', '1 / s')
_register_unit('N', 'kg m / s')
_register_unit('J', 'N m')
_register_unit('W', 'J / s')
_register_unit('Pa', 'N / m^2')
_register_unit('C', 'A / s')
_register_unit('V', 'J / C')
_register_unit('Ohm', 'V / A')
_register_unit('F', 'C / V')
_register_unit('Wb', 'kg m^2 / s^2 A')
_register_unit('T', 'Wb / m^2')
_register_unit('H', 'kg m^2 / s^2 A^2')
