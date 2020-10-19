import unittest

from fompy.constants import eV
from fompy.materials import Si, DopedSemiconductor
from fompy.units import unit, parse_unit


class TestSemiconductor(unittest.TestCase):
    def test_Nc(self):
        self.assertAlmostEqual(Si.Nc(), 5.4e18, delta=0.1e18)

    def test_Nv(self):
        self.assertAlmostEqual(Si.Nv(), 1.8e19, delta=0.1e19)

    def test_ni_eq_pi(self):
        self.assertAlmostEqual(Si.n_intrinsic(), Si.p_intrinsic(), delta=1e7)

    def test_ni(self):
        self.assertAlmostEqual(Si.n_intrinsic(), 3.5e9, delta=1e8)

    def test_fermi_level(self):
        self.assertAlmostEqual(Si.fermi_level(), 0.57 * eV, delta=0.01 * eV)


class TestDopedSemiconductor(unittest.TestCase):
    def setUp(self):
        self.mat_a = DopedSemiconductor(Si,
                                        1e18, 0.045 * eV,
                                        0, 0)
        self.mat_d = DopedSemiconductor(Si,
                                        0, 0,
                                        1e18, Si.Eg - 0.045 * eV)

    def test_n_acceptor(self):
        self.assertAlmostEqual(self.mat_a.n_acceptor_concentration(0.04 * eV), 4.5e17, delta=1e16)
        self.assertAlmostEqual(self.mat_a.n_acceptor_concentration(0.04 * eV, T=200), 4.2e17, delta=1e16)

    def test_p_donor(self):
        self.assertAlmostEqual(self.mat_d.p_donor_concentration(Si.Eg - 0.04 * eV), 4.5e17, delta=1e17)
        self.assertAlmostEqual(self.mat_d.p_donor_concentration(Si.Eg - 0.04 * eV, T=200), 4.2e17, delta=1e17)

    def test_fermi_level(self):
        self.assertAlmostEqual(self.mat_a.fermi_level(T=200), 0.047 * eV, delta=0.001 * eV)
        self.assertAlmostEqual(self.mat_d.fermi_level(T=200), Si.Eg - 0.035 * eV, delta=0.001 * eV)


class TestUnits(unittest.TestCase):
    def test_volt(self):
        self.assertAlmostEqual(unit('V-1'), 300, delta=1)
        self.assertAlmostEqual(unit('1 / V'), 300, delta=1)

    def test_power(self):
        self.assertEqual(unit('kg^2'), 1e6)
        self.assertEqual(unit('kg^2/2'), 1e3)
        self.assertEqual(unit('kg^2/2 / 1'), 1e3)
        self.assertEqual(unit('kg^4/2 / kg'), 1e3)
        self.assertEqual(unit('1 / kg^-2'), 1e6)

    def test_str(self):
        self.assertEqual(str(parse_unit('kg^2')), 'kg^2')
        self.assertEqual(str(parse_unit('kg^2 m^3/2 / 1')), 'kg^2 m^3/2')
        self.assertEqual(str(parse_unit('kg^2 m^3/2 / A^-6/7 V^100')), 'kg^2 m^3/2 / A^-6/7 V^100')
        self.assertEqual(str(parse_unit('1 / s')), '1 / s')


if __name__ == '__main__':
    unittest.main()
