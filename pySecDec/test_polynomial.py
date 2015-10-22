"""Unit tests for the Polynomial container class"""

from .polynomial import *
import unittest

class TestPolynomial(unittest.TestCase):
    def test_init(self):
        # Proper instantiation
        Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])

        # exponents should be integer
        self.assertRaisesRegexp(TypeError, "(e|E)xpolist.*integer", Polynomial, [(0,1.5),(1,0),(2,1)], ['A','B','C'])

        # Length mismatch between coeffs and expolist
        self.assertRaisesRegexp(AssertionError, "same length", Polynomial, [(0,1),(1,0),(2,1)], ['A','B','C','D'])


    def test_string_form(self):
        polynomial1 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'])
        string_polynomial1 = " + A*x0^0*x1^1 + B*x0^1*x1^0 + C*x0^2*x1^1 + D*x0^0*x1^0"
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)


    def test_copy(self):
        polynomial1 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'])
        polynomial2 = polynomial1.copy()

        self.assertEqual(str(polynomial1), str(polynomial2))

        polynomial1.expolist[0,0] = 5
        self.assertEqual(polynomial1.expolist[0,0],5)
        self.assertEqual(polynomial2.expolist[0,0],0)

        polynomial1.coeffs[0] += 'foo'
        self.assertEqual(polynomial1.coeffs[0],'Afoo')
        self.assertEqual(polynomial2.coeffs[0],'A')
