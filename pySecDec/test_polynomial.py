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
        string_polynomial1 = " + (A)*x1 + (B)*x0 + (C)*x0**2*x1 + (D)"
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)

        polynomial2 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'], polysymbols='z')
        string_polynomial2 = string_polynomial1.replace('x','z')
        self.assertEqual(str(polynomial2), string_polynomial2)
        self.assertEqual(repr(polynomial2), string_polynomial2)

        polynomial3 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'], polysymbols=['x','y'])
        string_polynomial3 = string_polynomial1.replace('x0','x').replace('x1','y')
        self.assertEqual(str(polynomial3), string_polynomial3)
        self.assertEqual(repr(polynomial3), string_polynomial3)

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
        self.assertEqual(str(polynomial), " + (2)*x0*x1**2 + (8)*x0**2*x1**2 + (2)*x0**2*x1 + (2)*x0**3*x1 + (6)*x0**3*x1**2")

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
        self.assertEqual(str(polysum), " + (1)*x1 + (7)*x0 + (10)*x0*x1 + (5)")

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
        self.assertEqual(str(p), ' + (a)*x0 + (b)*x1')
        self.assertEqual(str(p_squared), ' + (a**2)*x0**2 + (2*a*b)*x0*x1 + (b**2)*x1**2')

        p_sum = p_squared + SNCPolynomial([(1,1),(0,1)],[a,b])
        self.assertEqual(str(p_sum), ' + (a**2)*x0**2 + (2*a*b + a)*x0*x1 + (b**2)*x1**2 + (b)*x1')

    def test_empty_expolist(self):
        polynomial = SNCPolynomial([(0,1),(1,0),(2,1),(0,0)],[0,0,0,0])
        polynomial.combine()
        self.assertGreater(len(polynomial.expolist), 0)
        self.assertGreater(len(polynomial.coeffs), 0)

        np.testing.assert_array_equal(polynomial.expolist, [[0,0]])
        np.testing.assert_array_equal(polynomial.coeffs, [0])

    def test_creation_from_expression(self):
        from sympy import symbols, sympify
        x,y, a,b,c = symbols('x y a b c')
        polynomial_expression = a*x + b*y + c*x**2*y

        poly1 = SNCPolynomial.from_expression(polynomial_expression, [x,y])
        poly2 = SNCPolynomial.from_expression('a*x + b*y + c*x**2*y', ['x','y'])

        self.assertEqual((sympify(str(poly1)) - polynomial_expression).simplify(), 0)
        self.assertEqual((sympify(str(poly2)) - polynomial_expression).simplify(), 0)

        self.assertRaisesRegexp(TypeError, "\'x\*y\' is not.*symbol", SNCPolynomial.from_expression, 'a*x + b*y + c*x**2*y', [x,x*y])
        self.assertRaisesRegexp(TypeError, "polysymbols.*at least one.*symbol", SNCPolynomial.from_expression, 'a*x + b*y + c*x**2*y', [])

class TestExponentiatedPolynomial(unittest.TestCase):
    def test_init(self):
        ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3])
        ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3], exponent='A + eps')

    def test_string_form(self):
        polynomial1 = ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3], polysymbols='y')
        string_polynomial1 = ' + (x)*y0*y1**2 + (2)*y0 + (3)*y0**2*y1'
        # if exponent is one, do not show it
        self.assertEqual(str(polynomial1), string_polynomial1)
        self.assertEqual(repr(polynomial1), string_polynomial1)

        polynomial2 = ExponentiatedPolynomial([(1,2),(1,0),(2,1)], ['x',2,3], polysymbols='y', exponent='A + eps')
        string_polynomial2 = '( + (x)*y0*y1**2 + (2)*y0 + (3)*y0**2*y1)**(A + eps)'
        self.assertEqual(str(polynomial2), string_polynomial2)
        self.assertEqual(repr(polynomial2), string_polynomial2)

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
