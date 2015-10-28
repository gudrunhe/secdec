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

        # entries of expolist have variable length
        self.assertRaisesRegexp(AssertionError, "expolist.*same length", Polynomial, [(0,1,2),(1,0),(2,1)], ['A','B','C'])


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

    def test_has_constant_term(self):
        self.assertTrue(Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D']).has_constant_term())
        self.assertTrue(Polynomial([(0,1),(0,0),(1,0),(2,1),(4,0)],['A','B','C','D','']).has_constant_term())

        self.assertFalse(Polynomial([(0,1),(1,0),(2,1),(4,0)],['A','B','C','D']).has_constant_term())
        self.assertFalse(Polynomial([(0,1),(2,1),(4,0)],['A','B','D']).has_constant_term())

    def test_becomes_zero_for(self):
        self.assertTrue(Polynomial([(0,1,1,0),(2,1,0,5)],['A','B']).becomes_zero_for([1,0]))
        self.assertTrue(Polynomial([(0,1),(0,1),(1,5),(0,1),(4,4)],['A','B','C','D','']).becomes_zero_for([1]))

        self.assertFalse(Polynomial([(0,1),(1,0),(2,1),(4,0)],['A','B','C','D']).becomes_zero_for([1]))
        self.assertFalse(Polynomial([(0,1,0,1),(2,1,0,5),(0,0,3,5)],['A','B','C']).becomes_zero_for([1,0]))

class TestPolynomialProduct(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])
        p2 = Polynomial([(1,0,1),(2,1,0),(0,2,1)],['G','H','I'])

        # mismatch in number of parameters
        self.assertRaisesRegexp(TypeError, 'same number of variables.*all.*factors', PolynomialProduct,p0,p2)

        self.assertRaisesRegexp(AssertionError, 'at least one factor', PolynomialProduct)

        # Proper instantiation
        prod = PolynomialProduct(p0,p1)

        # made a copy?
        self.assertEqual(prod.factors[0].expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        prod.factors[0].expolist[0,0] = 5
        self.assertEqual(prod.factors[0].expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)
