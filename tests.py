import unittest
from math import sqrt, pi

import numpy as np

from fompy import functions
from fompy.constants import eV, volt, angstrom, amu, ampere
from fompy.materials import Si
from fompy.models import MSJunction, ContactType, DopedSemiconductor, Metal, PNJunction, \
    PNJunctionFullDepletion, PrimitiveCubicLattice, DiamondLikeLattice, FaceCenteredCubicLattice, \
    BodyCenteredCubicLattice, conductivity, concentration, PNJunctionNonDegenerate, hydrogen_like_energy, \
    hydrogen_like_radius
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

    def test_intrinsic_concentrations(self):
        mat = DopedSemiconductor(Si, 1e15, 0, 1e17, Si.Eg)  # random parameters
        self.assertAlmostEqual(Si.n_concentration(), Si.p_concentration(), delta=1e7)
        self.assertAlmostEqual(Si.n_concentration(), Si.i_concentration(), delta=1e7)
        self.assertAlmostEqual(mat.i_concentration(), Si.i_concentration(), delta=1e7)

    def test_ni(self):
        self.assertAlmostEqual(Si.n_concentration(), 3.5e9, delta=1e8)

    def test_fermi_level(self):
        self.assertAlmostEqual(Si.fermi_level(), 0.57 * eV, delta=0.01 * eV)

    def test_degenerate_fermi_level(self):
        mat = DopedSemiconductor(Si, 1e21, 0.045 * eV, 0, Si.Eg)
        self.assertAlmostEqual(mat.fermi_level() / eV, -0.04, delta=0.001)  # Note the negative value
        mat = DopedSemiconductor(Si, 10, 0, 1e20, Si.Eg - 0.045 * eV)
        self.assertAlmostEqual(mat.fermi_level() / eV, 1.142, delta=0.001)

    def test_intrinsic_fermi_level(self):
        mat = DopedSemiconductor(Si, 1e15, 0, 1e17, Si.Eg)  # random parameters
        self.assertAlmostEqual(mat.intrinsic_fermi_level(), 0.57 * eV, delta=0.01 * eV)
        self.assertEqual(mat.intrinsic_fermi_level(), Si.fermi_level())

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
        with self.assertRaises(ValueError):
            self.mat_a.conductivity_type(T=300, Ef=1 * eV)

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

    def test_prefixes(self):
        self.assertEqual(unit('dam'), 1e3)
        self.assertEqual(unit('dm'), 1e1)
        self.assertEqual(unit('dA'), ampere / 10)

    def test_error(self):
        with self.assertRaises(SyntaxError):
            parse_unit('kg^2 2')

        with self.assertRaises(SyntaxError):
            parse_unit('kg^2 / cm /')

        with self.assertRaises(SyntaxError):
            parse_unit('kg^2 / cm [')


class TestMSJunction(unittest.TestCase):
    def test_delta_phi(self):
        c = MSJunction(Metal(4.1 * eV), DopedSemiconductor(Si, 1e18, 0.045 * eV, 0, Si.Eg))
        self.assertAlmostEqual(c.delta_phi(300), -0.98 * volt, delta=0.1 * volt)

    def test_contact_type(self):
        with self.assertRaises(NotImplementedError):
            c = MSJunction(Metal(4.1 * eV), Si)
            c.contact_type()
        # Al -- p-Si
        c = MSJunction(Metal(4.1 * eV), DopedSemiconductor(Si, 1e17, 0.045 * eV, 0, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.INVERSION)
        # Pt -- n-Si
        c = MSJunction(Metal(5.2 * eV), DopedSemiconductor(Si, 0, 0, 1e18, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.INVERSION)
        # Cs -- n-Si
        c = MSJunction(Metal(2.14 * eV), DopedSemiconductor(Si, 0, 0, 1e18, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.AUGMENTATION)
        # Imaginary -- p-Si
        c = MSJunction(Metal(4.8 * eV), DopedSemiconductor(Si, 1e17, 0.045 * eV, 0, Si.Eg))
        self.assertEqual(c.contact_type(), ContactType.DEPLETION)

    def test_schottky_barrier(self):
        c = MSJunction(Metal(4.1 * eV), DopedSemiconductor(Si, 1e18, 0.045 * eV, 0, Si.Eg))
        self.assertAlmostEqual(c.schottky_barrier() / volt, 0.05, delta=0.001)

    def test_full_depletion_width(self):
        # The parameters are selected to comply with the condition delta_phi = 0.5 volt
        c = MSJunction(Metal(4.65 * eV), DopedSemiconductor(Si, 0, 0, 1e17, Si.Eg))
        self.assertAlmostEqual(c.delta_phi() / volt, 0.5, delta=0.001)
        self.assertAlmostEqual(c.full_depletion_width() / unit('nm'), 81, delta=1)

    def test_debye_length(self):
        c = MSJunction(Metal(4.1 * eV), DopedSemiconductor(Si, 0, 0, 1e18, Si.Eg))
        self.assertAlmostEqual(c.debye_length() / unit('nm'), 4.4, delta=0.1)


class TestFunctions(unittest.TestCase):
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
        self.assertEqual(0.0, np.max(np.abs(ys - functions.fd1(xs))))

    def test_issue_11(self):
        # Should not raise OverflowError
        m = DopedSemiconductor(Si, 1e12, 0, 0, Si.Eg)
        m.fermi_level(T=11)

        self.assertEqual(functions.fermi(Si.Eg, 0, 11), 0)


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


class TestPNJunctionNonDegenerate(unittest.TestCase):
    def test_np(self):
        n = p = 1e17
        pn = PNJunctionNonDegenerate(Si, n, 0, p, Si.Eg)
        self.assertAlmostEqual(pn.pn(0), 1.2e19, delta=1e18)
        self.assertAlmostEqual(pn.p_n(0), 123.3, delta=0.1)
        self.assertAlmostEqual(pn.n_p(0), 121.6, delta=0.1)

    def test_j0(self):
        n = p = 1e17
        pn = PNJunctionNonDegenerate(Si, n, 0, p, Si.Eg)
        d_n = 36
        d_p = 12
        l_n = l_p = 1e-2

        j0p = pn.j0_p(d_p, l_p)
        j0n = pn.j0_n(d_n, l_n)
        self.assertAlmostEqual((j0p + j0n) / unit('A / cm^2'), 9.4e-14, delta=0.1e-14)

    def test_current(self):
        n = p = 1e17
        pn = PNJunctionNonDegenerate(Si, n, 0, p, Si.Eg)
        d_n = 36
        d_p = 12
        l_n = l_p = 1e-2

        voltage = 0.8 * volt
        jp = pn.current_p(d_p, l_p, voltage)
        jn = pn.current_n(d_n, l_n, voltage)
        self.assertAlmostEqual((jp + jn) / unit('A / cm2'), 2.5, delta=0.1)

        voltage = -100 * volt
        jp = pn.current_p(d_p, l_p, voltage)
        jn = pn.current_n(d_n, l_n, voltage)
        self.assertAlmostEqual((jp + jn) / unit('A / cm2'), -9.4e-14, delta=0.1e-14)  # -J_0


class TestFormulae(unittest.TestCase):
    def setUp(self):
        self.mobility = 500 * unit('cm2 / V s')

    def test_conductivity(self):
        resistivity = 1 / conductivity(1e18, self.mobility)
        self.assertAlmostEqual(resistivity / unit('Ohm cm'), 0.012, delta=0.001)

    def test_concentration(self):
        self.assertAlmostEqual(concentration(0.01 * unit('Ohm cm'), self.mobility), 1.2e18, delta=1e17)
        self.assertAlmostEqual(concentration(10000 * unit('Ohm cm'), self.mobility), 1.2e12, delta=1e11)

    def test_hydrogen_like_energy(self):
        self.assertAlmostEqual(hydrogen_like_energy(Si.eps, Si.me) / eV, 0.035, delta=0.001)
        self.assertAlmostEqual(hydrogen_like_energy(Si.eps, Si.mh) / eV, 0.080, delta=0.001)

    def test_hydrogen_like_radius(self):
        self.assertAlmostEqual(hydrogen_like_radius(Si.eps, Si.me) / angstrom, 17.2, delta=0.1)
        self.assertAlmostEqual(hydrogen_like_radius(Si.eps, Si.mh) / angstrom, 7.6, delta=9.1)


if __name__ == '__main__':
    unittest.main()
