"""Unit tests for the Polynomial container class"""

from .decomposition import *
from .polynomial import Polynomial, PolynomialProduct
from .sector import Sector
from . import configure
import numpy as np
import unittest

def setUp():
    configure.powsymbol('^')
    configure.coeffs_in_parentheses(False)

class TestPrimaryDecomposition(unittest.TestCase):
    def setUp(self):
        # Feynman parameters are the ti
        # input is part of the 1Loop box

        # F = -s12*t1*t3 - s23*t0*t2
        self.F = Polynomial([(0,1,0,1),(1,0,1,0)],["-s12","-s23"])

        # U = t0 + t1 + t2 + t3
        self.U = Polynomial([(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)],["","","",""])

        self.initial_sector = Sector([self.F,self.U])

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

        self.assertEqual(str(F_primary[0]), " + -s12*x0^1*x1^0*x2^1 + -s23*x0^0*x1^1*x2^0")
        self.assertEqual(str(F_primary[1]), " + -s12*x0^0*x1^0*x2^1 + -s23*x0^1*x1^1*x2^0")
        self.assertEqual(str(F_primary[2]), " + -s12*x0^0*x1^1*x2^1 + -s23*x0^1*x1^0*x2^0")
        self.assertEqual(str(F_primary[3]), " + -s12*x0^0*x1^1*x2^0 + -s23*x0^1*x1^0*x2^1")


        np.testing.assert_array_equal(U_primary[0].expolist, np.array([(0,0,0),(1,0,0),(0,1,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[1].expolist, np.array([(1,0,0),(0,0,0),(0,1,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[2].expolist, np.array([(1,0,0),(0,1,0),(0,0,0),(0,0,1)]))
        np.testing.assert_array_equal(U_primary[3].expolist, np.array([(1,0,0),(0,1,0),(0,0,1),(0,0,0)]))

class TestIterativeDecomposition(unittest.TestCase):
    def setUp(self):
        # Feynman parameters are the ti
        # input is part of the 1Loop box

        # F = -s12*t1 - s23*t0*t2
        self.F = Polynomial([(0,1,0),(1,0,1)],["-s12","-s23"])

        # U = 1 + t0 + t1 + t2
        self.U = Polynomial([(0,0,0),(1,0,0),(0,1,0),(0,0,1)],["","","",""])

        # initialize an "empty" Jacobian in the sense that it is
        # the constant Polynomial with unit constant
        self.Jacobian = Polynomial([(0,0,0)],[""])

        self.sector = Sector([self.F,self.U])

    def test_refactorize(self):
        prod = PolynomialProduct(self.Jacobian, Polynomial([(1,1,0),(1,0,1)],["-s12","-s23"]))

        self.assertEqual(str(prod.factors[0]), ' + 1*x0^0*x1^0*x2^0')
        self.assertEqual(str(prod.factors[1]), ' + -s12*x0^1*x1^1*x2^0 + -s23*x0^1*x1^0*x2^1')

        copy0 = prod.copy()
        copy1 = prod.copy()
        copy2 = prod.copy()

        refactorize(copy0) # refactorize all parameters -> should find the factorization of parameter 0
        refactorize(copy1,0) # refactorize parameter 0 -> should find a factorization
        refactorize(copy2,1) # refactorize parameter 1 -> should NOT find the factorization of parameter 0

        self.assertEqual(str(copy0.factors[0]), str(copy1.factors[0]))
        self.assertEqual(str(copy0.factors[1]), str(copy1.factors[1]))

        self.assertEqual(str(copy1.factors[0]), ' + 1*x0^1*x1^0*x2^0')
        self.assertEqual(str(copy1.factors[1]), ' + -s12*x0^0*x1^1*x2^0 + -s23*x0^0*x1^0*x2^1')

        self.assertEqual(str(copy2.factors[0]), ' + 1*x0^0*x1^0*x2^0')
        self.assertEqual(str(copy2.factors[1]), ' + -s12*x0^1*x1^1*x2^0 + -s23*x0^1*x1^0*x2^1')

    def test_remap_parameters(self):
        remap_parameters([1,2], self.Jacobian, self.F, self.U) # modification in place
        self.assertEqual(str(self.Jacobian), " + 1*x0^0*x1^1*x2^0")
        self.assertEqual(str(self.F), " + -s12*x0^0*x1^1*x2^0 + -s23*x0^1*x1^1*x2^1")
        self.assertEqual(str(self.U), " + 1*x0^0*x1^0*x2^0 + 1*x0^1*x1^0*x2^0 + 1*x0^0*x1^1*x2^0 + 1*x0^0*x1^1*x2^1")

    def test_iteration_step(self):
        F = PolynomialProduct(self.Jacobian, self.F)
        U = PolynomialProduct(self.Jacobian, self.F)

        subsectors = iteration_step(self.sector)

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
        self.assertEqual(str(s0.Jacobian), " + 1*x0^1*x1^0*x2^0")
        self.assertEqual(str(s0_F.factors[0]), " + 1*x0^1*x1^0*x2^0")
        self.assertEqual(str(s0_F.factors[1]), " + -s12*x0^0*x1^1*x2^0 + -s23*x0^0*x1^0*x2^1")
        self.assertEqual(str(s0_U.factors[0]), " + 1*x0^0*x1^0*x2^0")
        self.assertEqual(str(s0_U.factors[1]), " + 1*x0^0*x1^0*x2^0 + 1*x0^1*x1^0*x2^0 + 1*x0^1*x1^1*x2^0 + 1*x0^0*x1^0*x2^1")

        # s1 should be the output of `remap_parameters([1,0], self.Jacobian, self.F, self.U)`
        # collected in a `Sector`
        self.assertEqual(str(s1.Jacobian), " + 1*x0^0*x1^1*x2^0")
        self.assertEqual(str(s1_F.factors[0]), " + 1*x0^0*x1^1*x2^0")
        self.assertEqual(str(s1_F.factors[1]), " + -s12*x0^0*x1^0*x2^0 + -s23*x0^1*x1^0*x2^1")
        self.assertEqual(str(s1_U.factors[0]), " + 1*x0^0*x1^0*x2^0")
        self.assertEqual(str(s1_U.factors[1]), " + 1*x0^0*x1^0*x2^0 + 1*x0^1*x1^1*x2^0 + 1*x0^0*x1^1*x2^0 + 1*x0^0*x1^0*x2^1")
