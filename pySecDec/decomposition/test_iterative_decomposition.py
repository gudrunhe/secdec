from .iterative import *
from .common import Sector
from ..algebra import Polynomial, ExponentiatedPolynomial, Product
from ..misc import sympify_expression
import numpy as np
import sympy as sp

import unittest
from nose.plugins.attrib import attr

class TestPrimaryDecomposition(unittest.TestCase):
    def setUp(self):
        # Feynman parameters are the ti
        # input is part of the 1Loop box

        # F = -s12*t1*t3 - s23*t0*t2
        self.F = Polynomial([(0,1,0,1),(1,0,1,0)],["-s12","-s23"])

        # U = t0 + t1 + t2 + t3
        self.U = Polynomial([(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)],[1,1,1,1])

        self.initial_sector = Sector([self.F,self.U])

    #@attr('active')
    def test_one_variable(self):
        U = Polynomial([[1]], [  1  ], ['x0'])
        F = Polynomial([[2]], ['msq'], ['x0'])
        initial_sector = Sector([U,F])
        primary_sectors = primary_decomposition(initial_sector)

        target_decomposed_U = sympify_expression(1)
        target_decomposed_F = sympify_expression('msq')

        self.assertEqual( len(primary_sectors) , 1 )

        self.assertEqual(  ( sympify_expression(primary_sectors[0].cast[0]) - target_decomposed_U ).simplify() , 0   )
        self.assertEqual(  ( sympify_expression(primary_sectors[0].cast[1]) - target_decomposed_F ).simplify() , 0   )

    #@attr('active')
    def test_primary_decomposition(self):
        primary_sectors = primary_decomposition(self.initial_sector)

        F_primary = [sector.cast[0].factors[1] for sector in primary_sectors]
        U_primary = [sector.cast[1].factors[1] for sector in primary_sectors]

        # 4 variables => 4 primary sectors
        self.assertEqual(len(F_primary),4)
        self.assertEqual(len(U_primary),4)


        np.testing.assert_array_equal(F_primary[0].expolist, np.array([(1,0,1),(0,1,0)]))
        np.testing.assert_array_equal(F_primary[1].expolist, np.array([(0,0,1),(1,1,0)]))
        np.testing.assert_array_equal(F_primary[2].expolist, np.array([(0,1,1),(1,0,0)]))
        np.testing.assert_array_equal(F_primary[3].expolist, np.array([(0,1,0),(1,0,1)]))

        self.assertEqual(str(F_primary[0]), " + (-s12)*x1*x3 + (-s23)*x2")
        self.assertEqual(str(F_primary[1]), " + (-s12)*x3 + (-s23)*x0*x2")
        self.assertEqual(str(F_primary[2]), " + (-s12)*x1*x3 + (-s23)*x0")
        self.assertEqual(str(F_primary[3]), " + (-s12)*x1 + (-s23)*x0*x2")


        np.testing.assert_array_equal(U_primary[0].expolist, np.array([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[1].expolist, np.array([(1,0,0),(0,0,0),(0,1,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[2].expolist, np.array([(1,0,0),(0,1,0),(0,0,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[3].expolist, np.array([(1,0,0),(0,1,0),(0,0,1),(0,0,0)]))

    #@attr('active')
    def test_select_indices(self):
        # do not consider the first and the last parameters
        indices = [1,2]
        primary_sectors = primary_decomposition(self.initial_sector, indices)

        F_primary = [sector.cast[0].factors[1] for sector in primary_sectors]
        U_primary = [sector.cast[1].factors[1] for sector in primary_sectors]

        # 2 integration variables => 2 primary sectors
        self.assertEqual(len(F_primary),2)
        self.assertEqual(len(U_primary),2)


        np.testing.assert_array_equal(F_primary[0].expolist, np.array([(0,0,1),(1,1,0)]))
        np.testing.assert_array_equal(F_primary[1].expolist, np.array([(0,1,1),(1,0,0)]))

        self.assertEqual(str(F_primary[0]), " + (-s12)*x3 + (-s23)*x0*x2")
        self.assertEqual(str(F_primary[1]), " + (-s12)*x1*x3 + (-s23)*x0")


        np.testing.assert_array_equal(U_primary[0].expolist, np.array([(1,0,0),(0,0,0),(0,1,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[1].expolist, np.array([(1,0,0),(0,1,0),(0,0,0),(0,0,1)]))

    #@attr('active')
    def test_keep_type(self):
        exponentiated_U = ExponentiatedPolynomial(self.U.expolist, self.U.coeffs, 'expoF', 'z')
        exponentiated_F = ExponentiatedPolynomial(self.F.expolist, self.F.coeffs, 'expoF', 'z')
        initial_sector = Sector([exponentiated_F,exponentiated_U])

        primary_sectors = primary_decomposition(initial_sector)

        F_primary = [sector.cast[0].factors[1] for sector in primary_sectors]
        U_primary = [sector.cast[1].factors[1] for sector in primary_sectors]

        # 4 variables => 4 primary sectors
        self.assertEqual(len(F_primary),4)
        self.assertEqual(len(U_primary),4)


        for poly in F_primary + U_primary:
            self.assertTrue(type(poly) is ExponentiatedPolynomial)

        np.testing.assert_array_equal(F_primary[0].expolist, np.array([(1,0,1),(0,1,0)]))
        np.testing.assert_array_equal(F_primary[1].expolist, np.array([(0,0,1),(1,1,0)]))
        np.testing.assert_array_equal(F_primary[2].expolist, np.array([(0,1,1),(1,0,0)]))
        np.testing.assert_array_equal(F_primary[3].expolist, np.array([(0,1,0),(1,0,1)]))

        self.assertEqual(str(F_primary[0]), "( + (-s12)*z1*z3 + (-s23)*z2)**(expoF)")
        self.assertEqual(str(F_primary[1]), "( + (-s12)*z3 + (-s23)*z0*z2)**(expoF)")
        self.assertEqual(str(F_primary[2]), "( + (-s12)*z1*z3 + (-s23)*z0)**(expoF)")
        self.assertEqual(str(F_primary[3]), "( + (-s12)*z1 + (-s23)*z0*z2)**(expoF)")


        np.testing.assert_array_equal(U_primary[0].expolist, np.array([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[1].expolist, np.array([(1,0,0),(0,0,0),(0,1,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[2].expolist, np.array([(1,0,0),(0,1,0),(0,0,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[3].expolist, np.array([(1,0,0),(0,1,0),(0,0,1),(0,0,0)]))

#@attr('active')
class TestIterativeDecomposition(unittest.TestCase):
    def setUp(self):
        # Feynman parameters are the ti
        # input is part of the 1Loop box

        # F = -s12*t1 - s23*t0*t2
        self.F = Polynomial([(0,1,0),(1,0,1)],["-s12","-s23"])

        # U = 1 + t0 + t1 + t2
        self.U = Polynomial([(0,0,0),(1,0,0),(0,1,0),(0,0,1)],[1,1,1,1])

        # initialize an "empty" Jacobian in the sense that it is
        # the constant Polynomial with unit constant
        self.Jacobian = Polynomial([(0,0,0)],[1])

        self.sector = Sector([self.F,self.U])

    def test_remap_parameters(self):
        remap_parameters([1,2], self.Jacobian, self.F, self.U) # modification in place
        self.assertEqual(str(self.Jacobian), " + (1)*x1")
        self.assertEqual(str(self.F), " + (-s12)*x1 + (-s23)*x0*x1*x2")
        self.assertEqual(str(self.U), " + (1) + (1)*x0 + (1)*x1 + (1)*x1*x2")

    #@attr('active')
    def test_find_singular_set(self):
        variables = ['x%i' %i for i in range(3)]

        F1 = Polynomial.from_expression('-s12*x1*x2**2 - s23*x0*x2 - s12*x0*x1',variables)
        sector1 = Sector([F1,self.U])
        self.assertEqual( find_singular_set(sector1),(0,2) )
        self.assertEqual( find_singular_set(sector1,indices=[0,1]),(0,1) )

        F2 = Polynomial.from_expression('-s12*x1**2*x2 - s23*x0*x2 - s12*x0*x1',variables)
        sector2 = Sector([F2,self.U])
        self.assertEqual( find_singular_set(sector2),(0,1) )

    def test_iteration_step(self):
        subsectors = list( iteration_step(self.sector) )

        # The algorithm should choose t0,t1
        # That generates two subsectors:
        self.assertEqual(len(subsectors), 2)

        s0 = subsectors[0]
        s0_F = s0.cast[0]
        s0_U = s0.cast[1]

        s1 = subsectors[1]
        s1_F = s1.cast[0]
        s1_U = s1.cast[1]

        # s0 should be the output of `remap_parameters([0,1], self.Jacobian, self.F, self.U)`
        # collected in a `Sector`
        self.assertEqual(str(s0.Jacobian), " + (1)*x0")
        self.assertEqual(str(s0_F.factors[0]), " + (1)*x0")
        self.assertEqual(str(s0_F.factors[1]), " + (-s12)*x1 + (-s23)*x2")
        self.assertEqual(str(s0_U.factors[0]), " + (1)")
        self.assertEqual(str(s0_U.factors[1]), " + (1) + (1)*x0 + (1)*x0*x1 + (1)*x2")

        # s1 should be the output of `remap_parameters([1,0], self.Jacobian, self.F, self.U)`
        # collected in a `Sector`
        self.assertEqual(str(s1.Jacobian), " + (1)*x1")
        self.assertEqual(str(s1_F.factors[0]), " + (1)*x1")
        self.assertEqual(str(s1_F.factors[1]), " + (-s12) + (-s23)*x0*x2")
        self.assertEqual(str(s1_U.factors[0]), " + (1)")
        self.assertEqual(str(s1_U.factors[1]), " + (1) + (1)*x0*x1 + (1)*x1 + (1)*x2")

    def test_iteration(self):
        subsectors = list( iterative_decomposition(self.sector) )

        # The algorithm should first choose t0,t1
        # That generates two subsectors; see test case above.
        # Next, it should select t1,t2 in sector 0 (see above) to further decompose F,
        # which increases the number of subsectors by one.
        # Sector 1 is already in standard form, it should not be further decomposed.
        self.assertEqual(len(subsectors), 3)

        s0_0 = subsectors[0]
        s0_0_F = s0_0.cast[0]
        s0_0_U = s0_0.cast[1]

        s0_1 = subsectors[1]
        s0_1_F = s0_1.cast[0]
        s0_1_U = s0_1.cast[1]

        s1 = subsectors[2]
        s1_F = s1.cast[0]
        s1_U = s1.cast[1]

        # sector 0_0
        self.assertEqual(str(s0_0.Jacobian), " + (1)*x0*x1")
        self.assertEqual(str(s0_0_F.factors[0]), " + (1)*x0*x1")
        self.assertEqual(str(s0_0_F.factors[1]), " + (-s12) + (-s23)*x2")
        self.assertEqual(str(s0_0_U.factors[0]), " + (1)")
        self.assertEqual(str(s0_0_U.factors[1]), " + (1) + (1)*x0 + (1)*x0*x1 + (1)*x1*x2")

        # sector 0_1
        self.assertEqual(str(s0_1.Jacobian), " + (1)*x0*x2")
        self.assertEqual(str(s0_1_F.factors[0]), " + (1)*x0*x2")
        self.assertEqual(str(s0_1_F.factors[1]), " + (-s12)*x1 + (-s23)")
        self.assertEqual(str(s0_1_U.factors[0]), " + (1)")
        self.assertEqual(str(s0_1_U.factors[1]), " + (1) + (1)*x0 + (1)*x0*x1*x2 + (1)*x2")

        # sector 1
        self.assertEqual(str(s1.Jacobian), " + (1)*x1")
        self.assertEqual(str(s1_F.factors[0]), " + (1)*x1")
        self.assertEqual(str(s1_F.factors[1]), " + (-s12) + (-s23)*x0*x2")
        self.assertEqual(str(s1_U.factors[0]), " + (1)")
        self.assertEqual(str(s1_U.factors[1]), " + (1) + (1)*x0*x1 + (1)*x1 + (1)*x2")

    #@attr('active')
    def test_multiple_iterative_decomposition(self):
        variables = ['x0', 'x1', 'x2']
        x0, x1, x2 = [Polynomial.from_expression('x%i' % i, variables) for i in range(3)]
        p0_mono = Polynomial.from_expression('x0*x1', variables)
        p0_poly = Polynomial([[1,0,0],[0,1,0],[0,0,1]],[1,1,1],variables) # "x0 + x1 + x2"
        p1 = Polynomial.from_expression('x1 + x2', variables)
        Jacobian = Polynomial.from_expression('x0', variables) # nonempty Jacobian to check transformation there
        initial_sector = Sector(cast=[Product(p0_mono,p0_poly)], other=[p1], Jacobian=Jacobian)
        psd_decomposition = list( iterative_decomposition(initial_sector) )

        target_decomposition = \
            [
                Sector(cast=[Product(p0_mono*x0**2,Polynomial([[0,0,0],[0,1,0],[0,0,1]],[1,1,1],variables))], other=[Polynomial([[1,1,0],[1,0,1]],[1,1], variables)], Jacobian=Jacobian*x0**2),
                Sector(cast=[Product(p0_mono*x1**2,Polynomial([[1,0,0],[0,0,0],[0,0,1]],[1,1,1],variables))], other=[Polynomial([[0,1,0],[0,1,1]],[1,1], variables)], Jacobian=Jacobian*x1*x1**2),
                Sector(cast=[Product(p0_mono*x2**3,Polynomial([[1,0,0],[0,1,0],[0,0,0]],[1,1,1],variables))], other=[Polynomial([[0,1,1],[0,0,1]],[1,1], variables)], Jacobian=Jacobian*x2*x2**2)
            ]

        self.assertEqual(len(psd_decomposition), 3)

        for stringify in (str, repr):
            self.assertEqual(stringify(psd_decomposition), stringify(target_decomposition))

    #@attr('active')
    def test_iterative_decomposition_selected_indices(self):
        variables = ['dummy0', 'x0', 'x1', 'dummy1', 'x2', 'eps']
        integration_variable_indices = [1,2,4]
        x0, x1, x2 = [Polynomial.from_expression('x%i' % i, variables) for i in range(3)]
        p0_mono = Polynomial.from_expression('x0*x1', variables)
        p0_poly = Polynomial([[0,1,0,0,0,0],[0,0,1,0,0,1],[0,0,0,0,1,0]],[1,1,1],variables) # "x0 + eps*x1 + x2"
        p1 = Polynomial.from_expression('x1 + eps*x2', variables)
        Jacobian = Polynomial.from_expression('eps*x0', variables)
        initial_sector = Sector(cast=[Product(p0_mono,p0_poly)], other=[p1], Jacobian=Jacobian)
        psd_decomposition = list( iterative_decomposition(initial_sector,integration_variable_indices) )

        target_decomposition = \
            [
                Sector(
                          cast=[Product(p0_mono*x0**2,Polynomial([[0,0,0,0,0,0],[0,0,1,0,0,1],[0,0,0,0,1,0]],[1,1,1],variables))],
                          other=[Polynomial([[0,1,1,0,0,0],[0,1,0,0,1,1]],[1,1], variables)], Jacobian=Jacobian*x0**2
                      ),
                Sector(
                          cast=[Product(p0_mono*x1**2,Polynomial([[0,1,0,0,0,0],[0,0,0,0,0,1],[0,0,0,0,1,0]],[1,1,1],variables))],
                          other=[Polynomial([[0,0,1,0,0,0],[0,0,1,0,1,1]],[1,1], variables)], Jacobian=Jacobian*x1*x1**2
                      ),
                Sector(
                          cast=[Product(p0_mono*x2**3,Polynomial([[0,1,0,0,0,0],[0,0,1,0,0,1],[0,0,0,0,0,0]],[1,1,1],variables))],
                          other=[Polynomial([[0,0,1,0,1,0],[0,0,0,0,1,1]],[1,1], variables)], Jacobian=Jacobian*x2*x2**2
                      )
            ]

        self.assertEqual(len(psd_decomposition), 3)

        for stringify in (str, repr):
            self.assertEqual(stringify(psd_decomposition), stringify(target_decomposition))
