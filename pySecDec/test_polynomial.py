"""Unit tests for the Polynomial container class"""

from .polynomial import *
from . import configure
import unittest

def setUp():
    configure.reset()
    configure.powsymbol('^')
    configure.coeffs_in_parentheses(False)

def tearDown():
    configure.reset()

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
        configure.powsymbol('^')
        polynomial1 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'])
        string_polynomial1 = " + A*x0^0*x1^1 + B*x0^1*x1^0 + C*x0^2*x1^1 + D*x0^0*x1^0"
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)

        configure.powsymbol('**')
        string_polynomial1 = string_polynomial1.replace('^','**')
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)

        configure.polysymbol('z')
        string_polynomial1 = string_polynomial1.replace('x','z')
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)

        setUp() # reset configuration

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

class TestSNCPolynomial(unittest.TestCase):
    def setUp(self):
        self.p0 = SNCPolynomial([(0,1),(1,0),(1,1)],[1,1,3])
        self.p1 = SNCPolynomial([(1,1),(1,1),(2,1)],[1,1,2])
        self.p2 = SNCPolynomial([(0,0),(1,0),(1,1)],[5,6,7])

    def test_init(self):
        self.assertRaisesRegexp(AssertionError, "coeffs.*one.dimensional", SNCPolynomial, [(0,3),(1,2)],[[1,4],[2,3]])
        SNCPolynomial([(1,2),(1,0),(2,1)], [1,2,3])

    def test_combine(self):
        polynomial = SNCPolynomial([[1,2],[1,2],[2,2]  ,  [2,1],[2,1],[3,1]  ,  [2,2],[2,2],[3,2]], [1,1,2  ,  1,1,2  ,  3,3,6])
        polynomial.combine()
        self.assertEqual(str(polynomial), " + 2*x0^1*x1^2 + 8*x0^2*x1^2 + 2*x0^2*x1^1 + 2*x0^3*x1^1 + 6*x0^3*x1^2")

    def test_mul(self):
        self.assertRaisesRegexp(TypeError, "unsupported operand type\(s\) for \*: \'SNCPolynomial\' and \'Polynomial\'", lambda: self.p0 * Polynomial([(0,0)],['']))
        self.assertRaisesRegexp(AssertionError, "Number of varibales must be equal for both factors in \*", lambda: self.p0 * SNCPolynomial([(0,0,3)],[4]))

        prod = self.p0 * self.p1

        #                                                    expolist    coefficient
        np.testing.assert_array_equal(prod.expolist, np.array([[1,2],   #     1
        #                                                      [1,2],   #     1
                                                               [2,2],   #     2
                                                               [2,1],   #     1
        #                                                      [2,1],   #     1
                                                               [3,1],   #     2
        #                                                      [2,2],   #     3
        #                                                      [2,2],   #     3
                                                               [3,2]])) #     6

        np.testing.assert_array_equal(prod.coeffs, [2,  8  ,  2,  2  ,      6])

    def test_add(self):
        self.assertRaisesRegexp(TypeError, "unsupported operand type\(s\) for \+: \'SNCPolynomial\' and \'Polynomial\'", lambda: self.p0 + Polynomial([(0,0)],['']))
        self.assertRaisesRegexp(AssertionError, "Number of varibales must be equal for both polynomials in \+", lambda: self.p0 + SNCPolynomial([(0,0,3)],[4]))

        polysum = self.p0 + self.p2
        self.assertEqual(str(polysum), " + 1*x0^0*x1^1 + 7*x0^1*x1^0 + 10*x0^1*x1^1 + 5*x0^0*x1^0")

    def test_negation(self):
        neg_p2 = - self.p2
        np.testing.assert_array_equal(neg_p2.coeffs, -self.p2.coeffs)
        np.testing.assert_array_equal(neg_p2.expolist, self.p2.expolist)

    def test_subtract(self):
        polysum = self.p0 + self.p2
        polysum -= self.p2
        self.assertEqual(str(polysum),str(self.p0))

    def test_sympy_binding(self):
        from sympy import symbols
        a,b = symbols('a b')

        p = SNCPolynomial([(1,0),(0,1)],[a,b])
        p_squared = p * p
        self.assertEqual(str(p), ' + a*x0^1*x1^0 + b*x0^0*x1^1')
        self.assertEqual(str(p_squared), ' + a^2*x0^2*x1^0 + 2*a*b*x0^1*x1^1 + b^2*x0^0*x1^2')

        configure.coeffs_in_parentheses(True)
        p_sum = p_squared + SNCPolynomial([(1,1),(0,1)],[a,b])
        self.assertEqual(str(p_sum), ' + (a^2)*x0^2*x1^0 + (2*a*b + a)*x0^1*x1^1 + (b^2)*x0^0*x1^2 + (b)*x0^0*x1^1')

        configure.coeffs_in_parentheses(False)

    def test_empty_expolist(self):
        polynomial = SNCPolynomial([(0,1),(1,0),(2,1),(0,0)],[0,0,0,0])
        polynomial.combine()
        self.assertGreater(len(polynomial.expolist), 0)
        self.assertGreater(len(polynomial.coeffs), 0)

        np.testing.assert_array_equal(polynomial.expolist, [[0,0]])
        np.testing.assert_array_equal(polynomial.coeffs, [0])

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

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])

        orig = PolynomialProduct(p0,p1)
        copy = orig.copy()

        self.assertEqual(orig.factors[0].expolist[0,0],0)
        self.assertEqual(copy.factors[0].expolist[0,0],0)
        orig.factors[0].expolist[0,0] = 5
        self.assertEqual(orig.factors[0].expolist[0,0],5)
        self.assertEqual(copy.factors[0].expolist[0,0],0)
