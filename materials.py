from constants import me, eV

class Semiconductor:
    def __init__(self, _me, _mh, _gap, _chi):
        """
        _me -- effective mass of electron
        _mh -- effective mass of hole
        _gap -- energy gap
        _chi -- electron affinity
        """
        self.me = _me
        self.mh = _mh
        self.Eg = _gap
        self.chi = _chi
# TODO: add classes for other material types, such as metals

# Values at 300K
# http://www.ioffe.ru/SVA/NSM/Semicond/Si/index.html
Si = Semiconductor(0.36*me, 0.81*me, 1.12*eV, 4.05*eV)
# TODO: add more materials
