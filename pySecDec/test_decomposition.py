"""Unit tests for the Polynomial container class"""

from .decomposition import *
from .polynomial import Polynomial
import numpy as np
import unittest

class TestPrimaryDecomposition(unittest.TestCase):
    def setUp(self):
        # Feynman parameters are the ti
        # input is part of the 1Loop box

        # F = -s12*t1*t3 - s23*t0*t2
        self.F = Polynomial([(0,1,0,1),(1,0,1,0)],["-s12","-s23"])

        # U = t0 + t1 + t2 + t3
        self.U = Polynomial([(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)],["","","",""])

    def test_primary_decomposition(self):
        F_primary = primary_decomposition(self.F)
        U_primary = primary_decomposition(self.U)


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

    def test_step(self):
        decompose_step([1,2], self.Jacobian, self.F, self.U) # modification in place
        self.assertEqual(str(self.Jacobian), " + 1*x0^0*x1^1*x2^0")
        self.assertEqual(str(self.F), " + -s12*x0^0*x1^1*x2^0 + -s23*x0^1*x1^1*x2^1")
        self.assertEqual(str(self.U), " + 1*x0^0*x1^0*x2^0 + 1*x0^1*x1^0*x2^0 + 1*x0^0*x1^1*x2^0 + 1*x0^0*x1^1*x2^1")

    def test_iteration(self):
        assert False
