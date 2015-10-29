"""Unit tests for the Sector container class"""

from .sector import *
from .polynomial import Polynomial, PolynomialProduct
import unittest

class TestSector(unittest.TestCase):
    def test_init(self):
        # Feynman parameters are the ti
        # input is part of the 1Loop box

        # F = -s12*t1 - s23*t0*t2
        F = Polynomial([(0,1,0),(1,0,1)],["-s12","-s23"])

        # U = 1 + t0 + t1 + t2
        U = Polynomial([(0,0,0),(1,0,0),(0,1,0),(0,0,1)],["","","",""])

        # "empty" Jacobian in the sense that it is
        # the constant Polynomial with unit constant
        Jacobian = Polynomial([(0,0,0)],[""])

        other_polynomial = Polynomial([(1,0,0,5),(0,1,0,2),(0,0,1,1)],["","",""])

        self.assertRaisesRegexp(AssertionError, 'Jacobian.*monomial', Sector, [F], Jacobian=U)
        self.assertRaisesRegexp(AssertionError, 'number of variables.*equal', Sector, [F], [other_polynomial])
        self.assertRaisesRegexp(AssertionError, '(f|F)irst factor.*monomial', Sector, [PolynomialProduct(F,U)])
        self.assertRaisesRegexp(AssertionError, 'two factors', Sector, [PolynomialProduct(F,U,Jacobian)])
        self.assertRaisesRegexp(AssertionError, 'at least one', Sector, [])
        Sector([PolynomialProduct(Jacobian,F)])

        sector = Sector([F])
        self.assertEqual(str(sector.Jacobian), str(Jacobian))

    def test_access(self):
        poly = Polynomial([(0,1,2),(1,0,5),(1,2,3),(9,4,2)],['','A','C','g'])
        sector = Sector([poly])

        self.assertEqual(sector.other,[])
        self.assertEqual(len(sector.cast),1)
        self.assertEqual(str(sector.cast[0].factors[1]),str(poly))
