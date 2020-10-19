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
    def __init__(self, name, power):
        prefix = None
        multiplier = 1
        for p, n in PREFIXES.items():
            if len(p) < len(name) and name.startswith(p):
                name = name[len(p):]
                prefix = p
                multiplier = n

        self.name = name
        self.prefix = prefix
        self.multiplier = multiplier
        self.value = UNITS[name]
        self.power = power

    def get_number(self):
        return (self.multiplier * self.value) ** self.power

    def __str__(self):
        res = (self.prefix or '') + self.name
        if self.power != 1:
            res += '^' + str(self.power)
        return res


class CompositeUnit:
    def __init__(self, nom, denom):
        self.nom = nom
        self.denom = denom

    def get_number(self):
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
    return parse_unit(text).get_number()


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
