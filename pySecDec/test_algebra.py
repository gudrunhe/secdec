"""Unit tests for the algebra module"""

from .algebra import *
import sympy as sp
import unittest
from nose.plugins.attrib import attr

class TestFunction(unittest.TestCase):
    def setUp(self):
        self.polysymbols = ['x0', 'x1', 'x2', 'x3']
        self.arg0 = Polynomial.from_expression('x0', self.polysymbols)
        self.arg1 = Polynomial.from_expression('x1', self.polysymbols)
        self.simplifyable_arg = Polynomial([[0,0,1,1],[0,0,1,1],[1,0,0,0]], [1,1,1], self.polysymbols)
        self.f = Function('f', self.arg0, self.arg1)
        self.g = Function('g', self.arg0, self.arg1, self.simplifyable_arg)

    #@attr('active')
    def test_init(self):
        # do not use variables defined in `setUp` here because of the last test in this function
        polysymbols = ['x0', 'x1', 'x2', 'x3']
        arg0 = Polynomial.from_expression('x0', polysymbols)
        arg1 = Polynomial.from_expression('x1', polysymbols)
        wrong_arg = Polynomial.from_expression('x4', polysymbols + ['x4'])
        self.assertRaisesRegexp(AssertionError, 'same.*number of variables.*all arguments', Function, 'f', arg0, arg1, wrong_arg)
        f = Function('f', arg0, arg1)

        self.assertEqual(f.number_of_variables, len(polysymbols))
        self.assertEqual( str(f) , 'f( + (1)*x0, + (1)*x1)' )

        # made a copy?
        arg0.expolist += 1
        self.assertEqual( str(f) , 'f( + (1)*x0, + (1)*x1)' )

    #@attr('active')
    def test_copy(self):
        f = self.f.copy()
        f.arguments[0].expolist[:,0] = 0
        f.arguments[0].expolist[:,1] = 1
        self.assertNotEqual( str(self.f) , str(f) )
        self.assertEqual( str(f) , 'f( + (1)*x1, + (1)*x1)' )

    #@attr('active')
    def test_simplify(self):
        unsimplified_g = self.g.copy()
        str_unsimplified_g = str(unsimplified_g)

        simplified_g = unsimplified_g.simplify()
        str_simplified_g = str(simplified_g)

        self.assertNotEqual(str_simplified_g, str_unsimplified_g)
        self.assertEqual( (sp.sympify(str_simplified_g) - sp.sympify(str_unsimplified_g)).simplify() , 0 )
        self.assertEqual( str(self.simplifyable_arg.simplify()) , str(simplified_g.arguments[-1]) )

    #@attr('active')
    def test_derive(self):
        # simple derivative
        df_d0 = self.f.derive(0)
        self.assertEqual( (sp.sympify(df_d0) - sp.sympify('df_d0(x0,x1)')).simplify() , 0 )

        # g( + (1)*x0, + (1)*x1, + (1)*x2*x3 + (1)*x2*x3 + (1)*x0)

        # derivatives with chain rule
        dg_d2 = self.g.derive(2)
        self.assertEqual( (sp.sympify(dg_d2) - sp.sympify('dg_d2(x0,x1,2*x2*x3+x0)*2*x3')).simplify() , 0 )

        dg_d0 = self.g.derive(0)
        self.assertEqual( (sp.sympify(dg_d0) - sp.sympify('dg_d0(x0,x1,2*x2*x3+x0) + dg_d2(x0,x1,2*x2*x3+x0)')).simplify() , 0 )

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
        from sympy import symbols
        A, Afoo = symbols('A Afoo')

        polynomial1 = Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'])
        polynomial2 = polynomial1.copy()

        self.assertEqual(str(polynomial1), str(polynomial2))

        polynomial1.expolist[0,0] = 5
        self.assertEqual(polynomial1.expolist[0,0],5)
        self.assertEqual(polynomial2.expolist[0,0],0)

        polynomial1.coeffs[0] = Afoo
        self.assertEqual(polynomial1.coeffs[0],Afoo)
        self.assertEqual(polynomial2.coeffs[0],A)

    def test_has_constant_term(self):
        self.assertTrue(Polynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D']).has_constant_term())
        self.assertTrue(Polynomial([(0,1),(0,0),(1,0),(2,1),(4,0)],['A','B','C','D',1]).has_constant_term())

        self.assertFalse(Polynomial([(0,1),(1,0),(2,1),(4,0)],['A','B','C','D']).has_constant_term())
        self.assertFalse(Polynomial([(0,1),(2,1),(4,0)],['A','B','D']).has_constant_term())

    def test_becomes_zero_for(self):
        self.assertTrue(Polynomial([(0,1,1,0),(2,1,0,5)],['A','B']).becomes_zero_for([1,0]))
        self.assertTrue(Polynomial([(0,1),(0,1),(1,5),(0,1),(4,4)],['A','B','C','D',1]).becomes_zero_for([1]))

        self.assertFalse(Polynomial([(0,1),(1,0),(2,1),(4,0)],['A','B','C','D']).becomes_zero_for([1]))
        self.assertFalse(Polynomial([(0,1,0,1),(2,1,0,5),(0,0,3,5)],['A','B','C']).becomes_zero_for([1,0]))

    def test_derive(self):
        from sympy import sympify
        polynomial = Polynomial([(2,1),(0,1)],['A', 'B'])

        derivative_0 = sympify( str(polynomial.derive(0)) )
        target_derivative_0 = sympify('2*A*x0*x1')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

        derivative_1 = sympify( str(polynomial.derive(1)) )
        target_derivative_1 = sympify('A*x0**2 + B')
        self.assertEqual( (derivative_1 - target_derivative_1).simplify() , 0 )

class TestFormerSNCPolynomial(unittest.TestCase):
    def setUp(self):
        self.p0 = Polynomial([(0,1),(1,0),(1,1)],[1,1,3])
        self.p1 = Polynomial([(1,1),(1,1),(2,1)],[1,1,2])
        self.p2 = Polynomial([(0,0),(1,0),(1,1)],[5,6,7])

    def test_init(self):
        self.assertRaisesRegexp(AssertionError, "coeffs.*one.dimensional", Polynomial, [(0,3),(1,2)],[[1,4],[2,3]])
        Polynomial([(1,2),(1,0),(2,1)], [1,2,3])

    def test_simplify(self):
        polynomial = Polynomial([[1,2],[1,2],[2,2]  ,  [2,1],[2,1],[3,1]  ,  [2,2],[2,2],[3,2]], [1,1,2  ,  1,1,2  ,  3,3,6])
        polynomial.simplify()
        self.assertEqual( (sp.sympify(str(polynomial)) - sp.sympify(" + (2)*x0*x1**2 + (8)*x0**2*x1**2 + (2)*x0**2*x1 + (2)*x0**3*x1 + (6)*x0**3*x1**2")).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(polynomial.coeffs), 5)
        self.assertEqual(polynomial.expolist.shape, (5,2))

    def test_mul(self):
        self.assertRaisesRegexp(AssertionError, "Number of varibales must be equal for both factors in \*", lambda: self.p0 * Polynomial([(0,0,3)],[4]))

        poly = Polynomial([(0,2),(1,0)],[1,'C'])
        intprod = 5 * poly
        self.assertEqual( (sp.sympify(str(intprod)) - sp.sympify(' + (5)*x1**2 + (5*C)*x0')).simplify() , 0)
        # should have minimal number of terms
        self.assertEqual(len(intprod.coeffs), 2)
        self.assertEqual(intprod.expolist.shape, (2,2))

        prod = self.p0 * self.p1

        #                              expolist    coefficient
        target_prod_expolist = np.array([[1,2],   #     1
        #                                [1,2],   #     1
                                         [2,2],   #     2
                                         [2,1],   #     1
        #                                [2,1],   #     1
                                         [3,1],   #     2
        #                                [2,2],   #     3
        #                                [2,2],   #     3
                                         [3,2]])  #     6

        target_prod_coeffs = np.array([2,  8  ,  2,  2  ,      6])


        # the order of the terms does not matter but for array comparison, it must be fixed
        sortkey_target = argsort_2D_array(target_prod_expolist)
        sorted_target_prod_expolist = target_prod_expolist[sortkey_target]
        sorted_target_prod_coeffs = target_prod_coeffs[sortkey_target]

        sortkey_prod = argsort_2D_array(prod.expolist)
        sorted_prod_expolist = prod.expolist[sortkey_prod]
        sorted_prod_coeffs = prod.coeffs[sortkey_prod]

        np.testing.assert_array_equal(sorted_prod_expolist, sorted_target_prod_expolist)
        np.testing.assert_array_equal(sorted_prod_coeffs, sorted_target_prod_coeffs)

    def test_add(self):
        self.assertRaisesRegexp(AssertionError, "Number of varibales must be equal for both polynomials in \+", lambda: self.p0 + Polynomial([(0,0,3)],[4]))

        self.assertEqual( str(10 + Polynomial([(1,0)],[1])) , ' + (10) + (1)*x0' )
        self.assertEqual( str(10 + Polynomial([(0,0),(1,0)],[2,1])) , ' + (12) + (1)*x0' )

        polysum = self.p0 + self.p2
        self.assertEqual( (sp.sympify(str(polysum)) - sp.sympify(" + (1)*x1 + (7)*x0 + (10)*x0*x1 + (5)")).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(polysum.coeffs), 4)
        self.assertEqual(polysum.expolist.shape, (4,2))


    def test_negation(self):
        neg_p2 = - self.p2
        np.testing.assert_array_equal(neg_p2.coeffs, -self.p2.coeffs)
        np.testing.assert_array_equal(neg_p2.expolist, self.p2.expolist)

    def test_subtract(self):
        polysum = self.p0 + self.p2
        polysum -= self.p2
        self.assertEqual(str(polysum),str(self.p0))

    #@attr('active')
    def test_pow(self):
        self.assertRaisesRegexp(TypeError, "unsupported", lambda: self.p0 ** 'a')
        self.assertRaisesRegexp(TypeError, "unsupported", lambda: self.p0 ** sp.sympify('a'))
        self.assertRaisesRegexp(ValueError, "must.*nonnegative", lambda: self.p0 ** -1)

        target_p0_zeroth_power = 0 * self.p0 + 1
        p0_zeroth_power = self.p0 ** 0
        self.assertEqual(sp.sympify(target_p0_zeroth_power - p0_zeroth_power), 0)

        target_p0_third_power = self.p0 * self.p0 * self.p0
        p0_third_power = self.p0 ** 3
        self.assertEqual(sp.sympify(target_p0_third_power - p0_third_power), 0)

        target_p1_sixth_power = self.p1 * self.p1 * self.p1 * self.p1 * self.p1 * self.p1
        p1_sixth_power = self.p1 ** 6
        self.assertEqual(sp.sympify(target_p1_sixth_power - p1_sixth_power), 0)

    def test_sympy_binding(self):
        from sympy import symbols
        a,b = symbols('a b')

        p = Polynomial([(1,0),(0,1)],[a,b])
        p_squared = p * p
        self.assertEqual(str(p), ' + (a)*x0 + (b)*x1')
        self.assertEqual( (sp.sympify(str(p_squared)) - sp.sympify(' + (a**2)*x0**2 + (2*a*b)*x0*x1 + (b**2)*x1**2')).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(p_squared.coeffs), 3)
        self.assertEqual(p_squared.expolist.shape, (3,2))

        p_sum = p_squared + Polynomial([(1,1),(0,1)],[a,b])
        self.assertEqual( (sp.sympify(str(p_sum)) - sp.sympify(' + (a**2)*x0**2 + (2*a*b + a)*x0*x1 + (b**2)*x1**2 + (b)*x1')).simplify() , 0)

        # should have minimal number of terms
        self.assertEqual(len(p_sum.coeffs), 4)
        self.assertEqual(p_sum.expolist.shape, (4,2))


    def test_empty_expolist(self):
        polynomial = Polynomial([(0,1),(1,0),(2,1),(0,0)],[0,0,0,0])
        polynomial.simplify()
        self.assertGreater(len(polynomial.expolist), 0)
        self.assertGreater(len(polynomial.coeffs), 0)

        np.testing.assert_array_equal(polynomial.expolist, [[0,0]])
        np.testing.assert_array_equal(polynomial.coeffs, [0])

    def test_creation_from_expression(self):
        from sympy import symbols, sympify
        x,y, a,b,c = symbols('x y a b c')
        polynomial_expression = a*x + b*y + c*x**2*y

        poly1 = Polynomial.from_expression(polynomial_expression, [x,y])
        poly2 = Polynomial.from_expression('a*x + b*y + c*x**2*y', ['x','y'])

        self.assertEqual((sympify(str(poly1)) - polynomial_expression).simplify(), 0)
        self.assertEqual((sympify(str(poly2)) - polynomial_expression).simplify(), 0)

        self.assertRaisesRegexp(TypeError, "\'x\*y\' is not.*symbol", Polynomial.from_expression, 'a*x + b*y + c*x**2*y', [x,x*y])
        self.assertRaisesRegexp(TypeError, "polysymbols.*at least one.*symbol", Polynomial.from_expression, 'a*x + b*y + c*x**2*y', [])

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

    def test_copy(self):
        polynomial1 = ExponentiatedPolynomial([(0,1),(1,0),(2,1),(0,0)],['A','B','C','D'],exponent='eps')
        polynomial2 = polynomial1.copy()

        self.assertEqual(str(polynomial1), str(polynomial2))
        self.assertEqual(polynomial1.exponent, polynomial2.exponent)

    def test_derive(self):
        from sympy import sympify, symbols
        A, B = symbols('A B')
        polynomial = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=sympify('a + b*eps'))

        derivative_0 = sympify( str(polynomial.derive(0)) )
        target_derivative_0 = sympify('(a + b*eps)*(A*x0**2*x1 + B)**(a + b*eps - 1) * (2*A*x0*x1)')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

        derivative_1 = sympify( str(polynomial.derive(1)) )
        target_derivative_1 = sympify('(a + b*eps)*(A*x0**2*x1 + B)**(a + b*eps - 1) * (A*x0**2)')
        self.assertEqual( (derivative_1 - target_derivative_1).simplify() , 0 )


        polynomial = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=Polynomial.from_expression('a + b*x0',['x0','x1']))
        derivative_0 = sympify( str(polynomial.derive(0)) )
        target_derivative_0 = sympify('(a + b*x0)*(A*x0**2*x1 + B)**(a + b*x0 - 1) * (2*A*x0*x1)   +   (A*x0**2*x1 + B)**(a + b*x0)*b*log(A*x0**2*x1 + B)')
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

    #@attr('active')
    def test_simplify(self):
        A, B = sp.symbols('A B')
        # <something>**0 = 1
        polynomial_to_power_zero_polynomial_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=Polynomial.from_expression('0',['x0','x1'])).simplify()
        polynomial_to_power_zero_sympy_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=sp.sympify('x-x')).simplify()
        polynomial_to_power_zero_numerical_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=0).simplify()

        for p in (polynomial_to_power_zero_polynomial_exponent, polynomial_to_power_zero_sympy_exponent, polynomial_to_power_zero_numerical_exponent):
            self.assertEqual(p.exponent, 1)
            np.testing.assert_array_equal(p.coeffs, [1])
            np.testing.assert_array_equal(p.expolist, [[0,0]])

        # <something>**1 = <something>
        polynomial_to_power_one_polynomial_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=Polynomial([[0,0,0],[0,0,0]], [2, -1])).simplify()
        polynomial_to_power_one_sympy_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=sp.sympify('1/2+1/2')).simplify()
        polynomial_to_power_one_numerical_exponent = ExponentiatedPolynomial([(2,1),(0,0)],[A, B],exponent=1).simplify()

        for p in (polynomial_to_power_one_polynomial_exponent, polynomial_to_power_one_sympy_exponent, polynomial_to_power_one_numerical_exponent):
            self.assertTrue(type(p) is Polynomial)
            np.testing.assert_array_equal(p.coeffs, [A,B])
            np.testing.assert_array_equal(p.expolist, [(2,1),(0,0)])

class TestProduct(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])
        p2 = Polynomial([(1,0,1),(2,1,0),(0,2,1)],['G','H','I'])

        # mismatch in number of parameters
        self.assertRaisesRegexp(TypeError, 'same number of variables.*all.*factors', Product,p0,p2)

        self.assertRaisesRegexp(AssertionError, 'at least one factor', Product)

        # Proper instantiation
        prod = Product(p0,p1)

        # made a copy?
        self.assertEqual(prod.factors[0].expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        prod.factors[0].expolist[0,0] = 5
        self.assertEqual(prod.factors[0].expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])

        orig = Product(p0,p1)
        copy = orig.copy()

        self.assertEqual(orig.factors[0].expolist[0,0],0)
        self.assertEqual(copy.factors[0].expolist[0,0],0)
        orig.factors[0].expolist[0,0] = 5
        self.assertEqual(orig.factors[0].expolist[0,0],5)
        self.assertEqual(copy.factors[0].expolist[0,0],0)

    #@attr('active')
    def test_derive(self):
        from sympy import sympify, symbols
        A, B = symbols('A B')
        polynomial1 = Polynomial([(2,1),(0,0)],[A, B])
        polynomial2 = Polynomial([(0,0),(1,0)],[1, 2])
        prod = Product(polynomial1, polynomial2)

        derivative_0 = prod.derive(0)
        derivative_0_1 = sympify( str(derivative_0.derive(1)) )
        target_derivative_0_1 = sympify('(A*x0**2*x1 + B) * (1 + 2*x0)')
        target_derivative_0_1 = sympify('(2*A*x0*x1) * (1 + 2*x0) + 2*(A*x0**2*x1 + B)')
        target_derivative_0_1 = sympify('(2*A*x0) * (1 + 2*x0) + 2*(A*x0**2)')
        self.assertEqual( (derivative_0_1 - target_derivative_0_1).simplify() , 0 )

    def test_string_form(self):
        p0 = ExponentiatedPolynomial([(0,1)],['A'],exponent='exponent')
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        prod = Product(p0,p1)
        string_prod = '(( + (A)*x1)**(exponent)) * ( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(prod), string_prod)
        self.assertEqual(repr(prod), string_prod)

#@attr('active')
class TestPow(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])
        p2 = Polynomial([(1,0,1),(2,1,0),(0,2,1)],['G','H','I'])

        # mismatch in number of parameters
        self.assertRaisesRegexp(TypeError, 'same number of variables.*base.*exponent', Pow,p0,p2)

        # Proper instantiation
        exp = Pow(p0,p1)

        # made a copy?
        self.assertEqual(exp.base.expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        exp.base.expolist[0,0] = 5
        self.assertEqual(exp.base.expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

        self.assertEqual(exp.exponent.expolist[0,0],8)
        self.assertEqual(p1.expolist[0,0],8)
        exp.exponent.expolist[0,0] = 5
        self.assertEqual(exp.exponent.expolist[0,0],5)
        self.assertEqual(p1.expolist[0,0],8)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])

        orig = Pow(p0,p1)
        copy = orig.copy()

        self.assertEqual(orig.base.expolist[0,0],0)
        self.assertEqual(copy.base.expolist[0,0],0)
        orig.base.expolist[0,0] = 5
        self.assertEqual(orig.base.expolist[0,0],5)
        self.assertEqual(copy.base.expolist[0,0],0)

        self.assertEqual(orig.exponent.expolist[0,0],8)
        self.assertEqual(copy.exponent.expolist[0,0],8)
        orig.exponent.expolist[0,0] = 5
        self.assertEqual(orig.exponent.expolist[0,0],5)
        self.assertEqual(copy.exponent.expolist[0,0],8)

    #@attr('active')
    def test_simplify(self):
        zero = Polynomial.from_expression('0', ['x0','x1','x2'])
        base = Polynomial([(0,0,1),(0,1,0),(1,0,0)], ['a','b','c'])

        one = Pow(base, zero).simplify()

        self.assertTrue(type(one) is Polynomial)
        np.testing.assert_array_equal(one.coeffs, [1])
        np.testing.assert_array_equal(one.expolist, [[0,0,0]])

        poly = Pow(base, one).simplify()
        self.assertTrue(type(poly) is Polynomial)
        np.testing.assert_array_equal(poly.coeffs, base.coeffs)
        np.testing.assert_array_equal(poly.expolist, base.expolist)

    #@attr('active')
    def test_derive(self):
        from sympy import sympify, symbols
        A, B = symbols('A B')
        polynomial1 = Polynomial.from_expression('A*x0 + B*x1', ['x0','x1'])
        polynomial2 = Polynomial.from_expression('x1', ['x0','x1'])
        exp = Pow(polynomial1, polynomial2)

        derivative_0 = exp.derive(0)
        target_derivative_0 = sympify('(A*x0 + B*x1) ** (x1 - 1) * A*x1')
        self.assertEqual( (sympify(derivative_0) - target_derivative_0).simplify() , 0 )

        derivative_0_1 = sympify(derivative_0.derive(1))
        target_derivative_0_1 = sympify('A*x1 *( (A*x0 + B*x1) ** (x1 - 1) * (log(A*x0 + B*x1) + (x1 - 1)*B/(A*x0 + B*x1)) ) + (A*x0 + B*x1) ** (x1 - 1) * A')
        self.assertEqual( (derivative_0_1 - target_derivative_0_1).simplify() , 0 )

        # `_Expression` in exponent
        exp = Pow(polynomial1, Sum(polynomial1, polynomial2))
        derivative_0 = exp.derive(0)
        target_derivative_0 = sympify('(A*x0 + B*x1)**(A*x0 + B*x1 + x1) * ( A*log(A*x0 + B*x1) + (A*x0 + B*x1 + x1)/(A*x0 + B*x1)*A )')
        self.assertEqual( (sympify(derivative_0) - target_derivative_0).simplify() , 0 )

    def test_string_form(self):
        p0 = ExponentiatedPolynomial([(0,1)],['A'],exponent='exponent')
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        exp = Pow(p0,p1)
        string_pow = '(( + (A)*x1)**(exponent)) ** ( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(exp), string_pow)
        self.assertEqual(repr(exp), string_pow)

#@attr('active')
class TestLog(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])

        # Proper instantiation
        ln = Log(p0)

        # made a copy?
        self.assertEqual(ln.arg.expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        ln.arg.expolist[0,0] = 5
        self.assertEqual(ln.arg.expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])

        orig = Log(p0)
        copy = orig.copy()

        self.assertEqual(orig.arg.expolist[0,0],0)
        self.assertEqual(copy.arg.expolist[0,0],0)
        orig.arg.expolist[0,0] = 5
        self.assertEqual(orig.arg.expolist[0,0],5)
        self.assertEqual(copy.arg.expolist[0,0],0)

    def test_simplify(self):
        one = Polynomial.from_expression(1, ['x0','x1','x2'])

        zero = Log(one).simplify()

        self.assertTrue(type(one) is Polynomial)
        self.assertEqual(sp.sympify(one), 1)
        np.testing.assert_array_equal(one.coeffs, [1])
        np.testing.assert_array_equal(one.expolist, [[0,0,0]])

    def test_derive(self):
        polynomial = Polynomial.from_expression('A*x0 + B*x1', ['x0','x1'])
        ln = Log(polynomial)

        derivative_0 = ln.derive(0).simplify()
        target_derivative_0 = sp.sympify('1/(A*x0 + B*x1)*A')
        self.assertEqual( (sp.sympify(derivative_0) - target_derivative_0).simplify() , 0 )

        derivative_0_1 = sp.sympify(derivative_0.derive(1))
        target_derivative_0_1 = sp.sympify('-A * (A*x0 + B*x1)**(-2) * B')
        self.assertEqual( (derivative_0_1 - target_derivative_0_1).simplify() , 0 )

    def test_string_form(self):
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        ln = Log(p1)
        string_ln = 'log( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(ln), string_ln)
        self.assertEqual(repr(ln), string_ln)

class TestSum(unittest.TestCase):
    def test_init(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])
        p2 = Polynomial([(1,0,1),(2,1,0),(0,2,1)],['G','H','I'])

        # mismatch in number of parameters
        self.assertRaisesRegexp(TypeError, 'same number of variables.*all.*summands', Sum,p0,p2)

        self.assertRaisesRegexp(AssertionError, 'at least one summand', Sum)

        # Proper instantiation
        psum = Sum(p0,p1)

        # made a copy?
        self.assertEqual(psum.summands[0].expolist[0,0],0)
        self.assertEqual(p0.expolist[0,0],0)
        psum.summands[0].expolist[0,0] = 5
        self.assertEqual(psum.summands[0].expolist[0,0],5)
        self.assertEqual(p0.expolist[0,0],0)

    def test_copy(self):
        p0 = Polynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        p1 = Polynomial([(8,1),(1,5),(2,1)],['D','E','F'])

        orig = Sum(p0,p1)
        copy = orig.copy()

        self.assertEqual(orig.summands[0].expolist[0,0],0)
        self.assertEqual(copy.summands[0].expolist[0,0],0)
        orig.summands[0].expolist[0,0] = 5
        self.assertEqual(orig.summands[0].expolist[0,0],5)
        self.assertEqual(copy.summands[0].expolist[0,0],0)

    def test_string_form(self):
        p0 = ExponentiatedPolynomial([(0,1)],['A'],exponent='exponent')
        p1 = Polynomial([(8,1),(1,5),(2,1)],['B','C','D'])
        sum = Sum(p0,p1)
        string_sum = '(( + (A)*x1)**(exponent)) + ( + (B)*x0**8*x1 + (C)*x0*x1**5 + (D)*x0**2*x1)'

        self.assertEqual(str(sum), string_sum)
        self.assertEqual(repr(sum), string_sum)

    def test_derive(self):
        from sympy import symbols, sympify
        A, B = symbols('A B')

        p0 = ExponentiatedPolynomial([(0,1)],[A])
        p1 = ExponentiatedPolynomial([(2,1)],[B])
        psum = Sum(p0,p1)

        derivative_0 = sympify( psum.derive(0) )
        target_derivative_0 = sympify( '2*B*x0*x1' )
        self.assertEqual( (derivative_0 - target_derivative_0).simplify() , 0 )

class TestLogOfPolynomial(unittest.TestCase):
    def test_string_form(self):
        p0 = LogOfPolynomial([(0,1),(1,0),(2,1)],['A','B','C'])
        str_p0 = 'log( + (A)*x1 + (B)*x0 + (C)*x0**2*x1)'

        self.assertEqual(str(p0),str_p0)
        self.assertEqual(repr(p0),str_p0)

    def test_construct_from_expression(self):
        p1 = LogOfPolynomial.from_expression('D*x0**8*x1 + E*x0*x1**5 + F*x0*x0*x1',['x0','x1'])
        str_p1 = 'log( + (D)*x0**8*x1 + (E)*x0*x1**5 + (F)*x0**2*x1)'
        sympy_p1 = sp.sympify(str_p1)

        self.assertEqual(str(p1),repr(p1))
        self.assertEqual( sp.sympify(repr(p1)) - sympy_p1 , 0 )

    def test_derive(self):
        expr = LogOfPolynomial([(2,1),(0,1)],['A', 'B'])

        derivative_0 = expr.derive(0)
        sympified_derivative_0 = sp.sympify( str(derivative_0) )
        target_derivative_0 = sp.sympify('1/(A*x0**2*x1 + B*x1) * 2*A*x0*x1')
        self.assertEqual( (sympified_derivative_0 - target_derivative_0).simplify() , 0 )
        self.assertEqual(type(derivative_0), Product)
        self.assertEqual(len(derivative_0.factors), 2)
        self.assertEqual(type(derivative_0.factors[0]), ExponentiatedPolynomial)

        derivative_1 = sp.sympify( str(expr.derive(1)) )
        target_derivative_1 = sp.sympify('1/(A*x0**2*x1 + B*x1) * (A*x0**2 + B)')
        self.assertEqual( (derivative_1 - target_derivative_1).simplify() , 0 )

    def test_simplify(self):
        # log(1) = 0
        expr = LogOfPolynomial([(0,0),(0,0)],[1, 0], polysymbols=['x','y'])
        simplified_expr = expr.simplify()

        self.assertTrue(type(simplified_expr) is Polynomial)
        np.testing.assert_array_equal(simplified_expr.expolist, [[0,0]])
        np.testing.assert_array_equal(simplified_expr.coeffs, [0])

#@attr('active')
class TestInsertion(unittest.TestCase):
    def test_insert_value_polynomial(self):
        poly = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        replaced_poly = replace(poly,index=1,value=sp.sympify('1/2'))
        self.assertEqual( sp.sympify(str(replaced_poly)) - sp.sympify('A + B*x0 + C/2 + D/4') , 0 )
        self.assertEqual(replaced_poly.number_of_variables, 2)

        removed_poly = replace(poly,index=1,value=sp.sympify('1/2'), remove=True)
        self.assertEqual( sp.sympify(str(replaced_poly)) - sp.sympify('A + B*x0 + C/2 + D/4') , 0 )
        self.assertEqual(removed_poly.number_of_variables, 1)
        self.assertEqual(removed_poly.expolist.shape[1], 1)

    def test_insert_value_polynomial_product(self):
        poly0 = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        poly1 = Polynomial([(0,0),(5,0)],['E',1])
        prod = Product(poly0,poly1)
        replaced_prod = replace(prod,index=1,value=0)
        self.assertEqual( (sp.sympify(str(replaced_prod)) - sp.sympify('(A + B*x0) * (E + x0**5)')).simplify() , 0 )
        self.assertEqual(replaced_prod.number_of_variables, 2)

        removed_prod = replace(prod, index=1, value=0, remove=True)
        self.assertEqual( (sp.sympify(str(removed_prod)) - sp.sympify('(A + B*x0) * (E + x0**5)')).simplify() , 0 )
        self.assertEqual(removed_prod.number_of_variables, 1)
        self.assertEqual(removed_prod.factors[0].expolist.shape[1], 1)
        self.assertEqual(removed_prod.factors[1].expolist.shape[1], 1)

    def test_insert_value_polynomial_sum(self):
        poly0 = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        poly1 = Polynomial([(0,0),(5,0)],['E',1])
        polysum = Sum(poly0,poly1)
        replaced_sum = replace(polysum,index=1,value=0)
        self.assertEqual( (sp.sympify(str(replaced_sum)) - sp.sympify('A + B*x0 + E + x0**5')).simplify() , 0 )
        self.assertEqual(replaced_sum.number_of_variables, 2)

        removed_sum = replace(polysum,index=1,value=0, remove=True)
        self.assertEqual( (sp.sympify(str(removed_sum)) - sp.sympify('(A + B*x0) + (E + x0**5)')).simplify() , 0 )
        self.assertEqual(removed_sum.number_of_variables, 1)
        self.assertEqual(removed_sum.summands[0].expolist.shape[1], 1)
        self.assertEqual(removed_sum.summands[1].expolist.shape[1], 1)

    def test_insert_in_exponent(self):
        exponent = Polynomial([(0,0),(5,0)],['E',1])
        poly1 = ExponentiatedPolynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'],exponent)
        replaced = replace(poly1,index=0,value=0)
        self.assertEqual( (sp.sympify(str(replaced)) - sp.sympify('(A + C*x1 + D*x1**2)**(E)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)
        self.assertEqual(replaced.expolist.shape[1], 2)
        self.assertEqual(replaced.exponent.number_of_variables, 2)
        self.assertEqual(replaced.exponent.expolist.shape[1], 2)

        removed = replace(poly1, index=0, value=2, remove=True)
        self.assertEqual( (sp.sympify(str(removed)) - sp.sympify('(A + 2*B + C*x1 + D*x1**2)**(E + 2**5)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)
        self.assertEqual(removed.exponent.number_of_variables, 1)
        self.assertEqual(removed.expolist.shape[1], 1)
        self.assertEqual(removed.exponent.expolist.shape[1], 1)

    def test_insert_in_pow(self):
        exponent = Polynomial([(0,0),(5,0)],['E',1])
        base = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        expr = Pow(base, exponent)

        replaced = replace(expr,index=0,value=0)
        self.assertEqual( (sp.sympify(str(replaced)) - sp.sympify('(A + C*x1 + D*x1**2)**(E)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)

        removed = replace(expr, index=0, value=2, remove=True)
        self.assertEqual( (sp.sympify(str(removed)) - sp.sympify('(A + 2*B + C*x1 + D*x1**2)**(E + 2**5)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)

    def test_insert_in_log(self):
        arg = Polynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'])
        expr = Log(arg)

        replaced = replace(expr,index=1,value=sp.sympify('(3*k+1)'))
        self.assertEqual( (sp.sympify(str(replaced)) - sp.sympify('log(A + B*x0 + C*(3*k+1) + D*(3*k+1)**2)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)

        removed = replace(expr, index=1, value=sp.sympify('(3*k+1)'), remove=True)
        self.assertEqual( (sp.sympify(str(removed)) - sp.sympify('log(A + B*x0 + C*(3*k+1) + D*(3*k+1)**2)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)

    #@attr('active')
    def test_insert_in_Function(self):
        exponent = Polynomial([(0,0),(5,0)],['E',1])
        poly1 = ExponentiatedPolynomial([(0,0),(1,0),(0,1),(0,2)],['A','B','C','D'],exponent)
        poly2 = Polynomial([(1,0),(0,1)],[1,1])
        f = Function('f', poly1, poly2)

        replaced = replace(f,index=0,value=0)
        self.assertEqual( (sp.sympify(str(replaced)) - sp.sympify('f( (A + C*x1 + D*x1**2)**(E), x1)')).simplify() , 0 )
        self.assertEqual(replaced.number_of_variables, 2)
        self.assertEqual(replaced.arguments[0].expolist.shape[1], 2)
        self.assertEqual(replaced.arguments[1].expolist.shape[1], 2)
        self.assertEqual(replaced.arguments[1].number_of_variables, 2)
        self.assertEqual(replaced.arguments[0].exponent.number_of_variables, 2)
        self.assertEqual(replaced.arguments[0].exponent.expolist.shape[1], 2)

        removed = replace(f, index=0, value=2, remove=True)
        self.assertEqual( (sp.sympify(str(removed)) - sp.sympify('f( (A + 2*B + C*x1 + D*x1**2)**(E + 2**5), 2 + x1)')).simplify() , 0 )
        self.assertEqual(removed.number_of_variables, 1)
        self.assertEqual(removed.arguments[0].expolist.shape[1], 1)
        self.assertEqual(removed.arguments[1].expolist.shape[1], 1)
        self.assertEqual(removed.arguments[1].number_of_variables, 1)
        self.assertEqual(removed.arguments[0].exponent.number_of_variables, 1)
        self.assertEqual(removed.arguments[0].exponent.expolist.shape[1], 1)

    def test_error(self):
        self.assertRaisesRegexp(TypeError, 'Can.*only.*Polynomial.*not.*int', replace, 3, 2, 1)

#@attr('active')
class TestGetSymbols(unittest.TestCase):
    def setUp(self):
        self.polysymbols = sp.sympify(['x','y'])
        self.p0 = Polynomial.from_expression('a*x + b*y', self.polysymbols)
        self.p1 = Polynomial.from_expression('c*x + d*y', self.polysymbols)

    def test_get_from_polynomial(self):
        self.assertEqual(self.p0.symbols, self.polysymbols)
        self.assertEqual(self.p1.symbols, self.polysymbols)

    def test_get_from_sum(self):
        expr = Sum(self.p0, self.p1)
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_product(self):
        expr = Product(self.p0, self.p1)
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_log(self):
        expr = Sum(Log(self.p0), Log(self.p1))
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_pow(self):
        expr = Pow(Log(self.p0), Log(self.p1))
        self.assertEqual(expr.symbols, self.polysymbols)

    def test_get_from_function(self):
        expr = Function('f', self.p0, self.p1)
        self.assertEqual(expr.symbols, self.polysymbols)
