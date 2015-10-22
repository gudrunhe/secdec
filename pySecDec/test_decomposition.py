"""Unit tests for the Polynomial container class"""

from .decomposition import *
from .polynomial import Polynomial
import unittest

class TestPrimaryDecomposition(unittest.TestCase):
    def test_primary_decomposition(self):

        
        assert False

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
