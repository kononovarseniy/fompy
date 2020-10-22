import unittest
from math import sqrt, pi

import numpy as np

from fompy.constants import eV, volt, angstrom, amu
from fompy.functions import fd1
from fompy.materials import Si
from fompy.models import MetalSemiconductorContact, ContactType, DopedSemiconductor, Metal, PNJunction, \
    PNJunctionFullDepletion, PrimitiveCubicLattice, DiamondLikeLattice, FaceCenteredCubicLattice, \
    BodyCenteredCubicLattice
from fompy.units import unit, parse_unit


class TestCrystalLattice(unittest.TestCase):
    def test_getters(self):
        lat = PrimitiveCubicLattice(5.43 * angstrom, 28 * amu)
        self.assertEqual(lat.a, 5.43 * angstrom)
        self.assertEqual(lat.m, 28 * amu)

    def test_primitive(self):
        lat = PrimitiveCubicLattice(angstrom, amu)
        self.assertEqual(lat.r, angstrom / 2)
        self.assertEqual(lat.N, 1)

    def test_face_centered(self):
        lat = FaceCenteredCubicLattice(angstrom, amu)
        self.assertEqual(lat.r, angstrom * sqrt(2) / 4)
        self.assertEqual(lat.N, 4)

    def test_body_centered(self):
        lat = BodyCenteredCubicLattice(angstrom, amu)
        self.assertEqual(lat.r, angstrom * sqrt(3) / 4)
        self.assertEqual(lat.N, 2)

    def test_diamond_like(self):
        lat = DiamondLikeLattice(5.43 * angstrom, 28 * amu)
        self.assertEqual(lat.r, 5.43 * angstrom * sqrt(3) / 8)
        self.assertEqual(lat.N, 8)

        self.assertAlmostEqual(lat.concentration, 4.997e22, delta=0.001e22)
        self.assertAlmostEqual(lat.density, 2.324, delta=0.001)

    def test_packing_density(self):
        lp = PrimitiveCubicLattice(angstrom, amu)
        lf = FaceCenteredCubicLattice(angstrom, amu)
        lb = BodyCenteredCubicLattice(angstrom, amu)
        ld = DiamondLikeLattice(angstrom, amu)
        self.assertAlmostEqual(lp.packing_density, pi / 6)
        self.assertAlmostEqual(lf.packing_density, pi * sqrt(2) / 6)
        self.assertAlmostEqual(lb.packing_density, pi * sqrt(3) / 8)
        self.assertAlmostEqual(ld.packing_density, pi * sqrt(3) / 16)


class TestSemiconductor(unittest.TestCase):
    def test_Nc(self):
        self.assertAlmostEqual(Si.Nc(), 5.4e18, delta=0.1e18)

    def test_Nv(self):
        self.assertAlmostEqual(Si.Nv(), 1.8e19, delta=0.1e19)

    def test_ni_eq_pi(self):
        self.assertAlmostEqual(Si.n_concentration(), Si.p_concentration(), delta=1e7)

    def test_ni(self):
        self.assertAlmostEqual(Si.n_concentration(), 3.5e9, delta=1e8)

    def test_fermi_level(self):
        self.assertAlmostEqual(Si.fermi_level(), 0.57 * eV, delta=0.01 * eV)

    def test_conductivity_type(self):
        self.assertEqual(Si.conductivity_type(), 'i')


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

    def test_conductivity_type(self):
        self.assertEqual(self.mat_a.conductivity_type(), 'p')
        self.assertEqual(self.mat_d.conductivity_type(), 'n')


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


class TestMetalSemiconductorContact(unittest.TestCase):
    def test_delta_phi(self):
        c = MetalSemiconductorContact(Metal(4.1 * eV), DopedSemiconductor(Si, 1e18, 0.045 * eV, 0, Si.Eg))
        self.assertAlmostEqual(c.delta_phi(300), -0.98 * volt, delta=0.1 * volt)

    def test_contact_type(self):
        # Al -- p-Si
        c = MetalSemiconductorContact(Metal(4.1 * eV), DopedSemiconductor(Si, 1e17, 0.045 * eV, 0, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.INVERSION)
        # Pt -- n-Si
        c = MetalSemiconductorContact(Metal(5.2 * eV), DopedSemiconductor(Si, 0, 0, 1e18, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.INVERSION)
        # Cs -- n-Si
        c = MetalSemiconductorContact(Metal(2.14 * eV), DopedSemiconductor(Si, 0, 0, 1e18, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.AUGMENTATION)
        # Imaginary -- p-Si
        c = MetalSemiconductorContact(Metal(4.8 * eV), DopedSemiconductor(Si, 1e17, 0.045 * eV, 0, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.DEPLETION)


class TestFermiDiracIntegral(unittest.TestCase):
    def test_fd1(self):
        xs = np.array([-3.0, -1.0, 1.0, 3.0, 7.0, 15.0, 30.0, 60.0])
        ys = np.array([0.04336636755041557,
                       0.2905008961699176,
                       1.396375280666564,
                       3.976985354047977,
                       12.664637572252104,
                       38.9430466009327,
                       109.69481833726653,
                       309.9448732700438
                       ])
        self.assertEqual(0.0, np.max(np.abs(ys - fd1(xs))))


class TestPNJunction(unittest.TestCase):
    def setUp(self):
        self.pn = PNJunction(Si, 5e16, None, 1e16, None)
        self.pn_fd = PNJunctionFullDepletion(Si, 5e16, None, 1e16, None)
        self.pn_fd2 = PNJunctionFullDepletion(Si, 1e16, None, 1e16, None)

    def test_delta_phi(self):
        self.assertAlmostEqual(self.pn.delta_phi() / volt, 0.81, delta=0.05)

    def test_delta_phi_n(self):
        self.assertAlmostEqual(self.pn_fd.delta_phi_n() / volt, 0.71, delta=0.05)

    def test_delta_phi_p(self):
        self.assertAlmostEqual(self.pn_fd.delta_phi_p() / volt, 0.14, delta=0.05)

    def test_delta_w(self):
        self.assertAlmostEqual(self.pn_fd.w(), 3.5e-5, delta=1e-6)
        self.assertAlmostEqual(self.pn_fd2.w(), 4.4e-5, delta=1e-6)

    def test_delta_w_n(self):
        self.assertAlmostEqual(self.pn_fd.w_n(), 2.9e-5, delta=1e-6)
        self.assertAlmostEqual(self.pn_fd2.w_n(), 2.2e-5, delta=1e-6)

    def test_delta_w_p(self):
        self.assertAlmostEqual(self.pn_fd.w_p(), 5.9e-6, delta=1e-7)
        self.assertAlmostEqual(self.pn_fd2.w_p(), 2.2e-5, delta=1e-6)

    def test_neutrality(self):
        self.assertEqual(self.pn_fd.w_p() * self.pn_fd.p_mat.n_acceptor_concentration() -
                         self.pn_fd.w_n() * self.pn_fd.n_mat.p_donor_concentration(), 0.0)
        self.assertEqual(self.pn_fd2.w_p() * self.pn_fd2.p_mat.n_acceptor_concentration() -
                         self.pn_fd2.w_n() * self.pn_fd2.n_mat.p_donor_concentration(), 0.0)


if __name__ == '__main__':
    unittest.main()
