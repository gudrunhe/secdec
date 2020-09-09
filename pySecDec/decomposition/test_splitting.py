from .splitting import *
from .common import Sector
from ..algebra import Polynomial, ExponentiatedPolynomial
from ..misc import sympify_expression
import numpy as np
import sympy as sp
import unittest
from nose.plugins.attrib import attr

class TestRemapOneToZero(unittest.TestCase):
    def setUp(self):
        # U and F generated from:
        # propagators = {k1^2, (k1+p1+p2)^2, (k1-k2)^2, (k1-k2+p1)^2-s, (k2)^2, (k2+p2)^2};
        # Dim = 4-2*eps;
        # ExternalMomenta = {p1,p2,p3};
        # Masses = KinematicInvariants = {s};
        # ScalarProductRules = {
        #                           SP[p1,p1]->0,
        #                           SP[p2,p2]->0,
        #                           SP[p3,p3]->s,
        #                           SP[p1,p2]->s/2,
        #                           SP[p2,p3]->-s/2,
        #                           SP[p1,p3]->-s/2
        #                      };

        x0, x1, x2, x3, x4, x5 = self.Feynman_parameters = sp.symbols('x0, x1, x2, x3, x4, x5')
        s = self.s = sp.symbols('s')
        self.sympy_U =  + (1)*x3*x5 + (1)*x3*x4 + (1)*x2*x5 + (1)*x2*x4 + (1)*x1*x5 + (1)*x1*x4 \
                        + (1)*x1*x3 + (1)*x1*x2 + (1)*x0*x5 + (1)*x0*x4 + (1)*x0*x3 + (1)*x0*x2
        self.U = Polynomial.from_expression(self.sympy_U, self.Feynman_parameters)

        self.sympy_F = + (s)*x3**2*x5 + (s)*x3**2*x4 + (s)*x2*x3*x5 + (s)*x2*x3*x4 \
                       + (s)*x1*x3*x5 + (s)*x1*x3*x4 + (s)*x1*x3**2 + (-s)*x1*x2*x4 \
                       + (s)*x1*x2*x3 + (s)*x0*x3*x4 + (s)*x0*x3**2 + (s)*x0*x2*x3 \
                       + (-s)*x0*x1*x5 + (-s)*x0*x1*x4 + (-s)*x0*x1*x3 + (-s)*x0*x1*x2
        self.F = Polynomial.from_expression(self.sympy_F, self.Feynman_parameters)

        self.initial_sector = Sector([self.U,self.F])

    #@attr('active')
    def test_remap_one_to_zero_U(self):
        x0, x1, x2, x3, x4, x5 = self.Feynman_parameters
        s = self.s
        remap_indices = [1,5]
        remapped_U = sympify_expression( remap_one_to_zero(self.U, *remap_indices) )
        target_remapped_U = + (1)*x3*(1-x5) + (1)*x3*x4 + (1)*x2*(1-x5) + (1)*x2*x4 \
                            + (1)*(1-x1)*(1-x5) + (1)*(1-x1)*x4 + (1)*(1-x1)*x3 \
                            + (1)*(1-x1)*x2 + (1)*x0*(1-x5) + (1)*x0*x4 + (1)*x0*x3 \
                            + (1)*x0*x2
        self.assertEqual( (remapped_U - target_remapped_U).simplify() , 0 )

    #@attr('active')
    def test_remap_one_to_zero_F(self):
        x0, x1, x2, x3, x4, x5 = self.Feynman_parameters
        s = self.s
        remap_indices = [0,1,2]
        remapped_F = sympify_expression( remap_one_to_zero(self.F, *remap_indices) )
        target_remapped_F = + (s)*x3**2*x5 + (s)*x3**2*x4 + (s)*(1-x2)*x3*x5 + (s)*(1-x2)*x3*x4 \
                            + (s)*(1-x1)*x3*x5 + (s)*(1-x1)*x3*x4 + (s)*(1-x1)*x3**2 \
                            + (-s)*(1-x1)*(1-x2)*x4 + (s)*(1-x1)*(1-x2)*x3 + (s)*(1-x0)*x3*x4 \
                            + (s)*(1-x0)*x3**2 + (s)*(1-x0)*(1-x2)*x3 + (-s)*(1-x0)*(1-x1)*x5 \
                            + (-s)*(1-x0)*(1-x1)*x4 + (-s)*(1-x0)*(1-x1)*x3 + (-s)*(1-x0)*(1-x1)*(1-x2)
        self.assertEqual( (remapped_F - target_remapped_F).simplify() , 0 )

    #@attr('active')
    def test_remap_exponentiated(self):
        x0, x1, x2, x3, x4, x5 = self.Feynman_parameters
        polynomial = Polynomial.from_expression(x0 - x1, self.Feynman_parameters)
        exponent = sp.symbols('exponent')
        exponentiated_polynomial = ExponentiatedPolynomial(polynomial.expolist, polynomial.coeffs,
                                                           exponent, polynomial.polysymbols)

        remapped_exponentiated_polynomial = sympify_expression( remap_one_to_zero(exponentiated_polynomial, 0,1) )
        target_remapped_exponentiated_polynomial = (  (1-x0) - (1-x1)  ) ** exponent

        self.assertEqual(  (remapped_exponentiated_polynomial - target_remapped_exponentiated_polynomial).simplify() , 0  )

    #@attr('active')
    def test_remap_constant(self):
        polynomial = Polynomial.from_expression('coeff', self.Feynman_parameters)
        remapped = remap_one_to_zero(polynomial, 0,1)
        target_remapped = sympify_expression('coeff')

        self.assertTrue( type(remapped) is Polynomial )
        self.assertEqual(  (sympify_expression(remapped) - target_remapped).simplify() , 0  )

class TestFindSingularSetsAtOne(unittest.TestCase):
    #@attr('active')
    def test_find_singular_sets_at_one_empty(self):
        poly = Polynomial.from_expression('x0 - x1', ['x0','x1'])
        singular_set = find_singular_sets_at_one(poly)
        self.assertEqual(singular_set, [tuple(),(0,1)])

    #@attr('active')
    def test_find_singular_sets_at_one_simple(self):
        poly = Polynomial.from_expression('1 - x0 + x1', ['x0','x1'])
        singular_set = find_singular_sets_at_one(poly)
        self.assertEqual(singular_set, [(0,)])

    #@attr('active')
    def test_find_singular_sets_at_one_medium(self):
        poly = Polynomial.from_expression('2 - x0 + a*x1**2*x5**4 - x3*x0**8', ['x0','x1','x2','x3','x4','x5'])
        singular_set = find_singular_sets_at_one(poly)
        self.assertEqual(singular_set, [(0,3),(0,1,3),(0,3,5)])

    #@attr('active')
    def test_find_singular_sets_at_one_complicated(self):
        s, x0, x1, x2, x3, x4, x5 = sp.symbols('s, x0, x1, x2, x3, x4, x5')
        poly = + (-s)*x0*x2 + (-s)*x0*x1*x3 + (-s)*x4 + (-s) + (s)*x0*x2*x3 \
               + (s)*x0*x1*x3**2 + (s)*x3*x4 + (s)*x0*x1*x2*x3 + (-s)*x2*x4 \
               + (s)*x0*x1**2*x3**2 + (s)*x1*x3*x4 + (s)*x1*x3 + (s)*x2*x3*x4 \
               + (s)*x2*x3 + (s)*x1*x3**2*x4 + (s)*x1*x3**2
        poly = Polynomial.from_expression(poly, ['x0','x1','x2','x3','x4','x5'])
        singular_set = find_singular_sets_at_one(poly)
        self.assertEqual(singular_set, [(2,3),(0,2,3),(2,3,4),(0,2,3,4)])

    #@attr('active')
    def test_find_singular_sets_only_at_one(self):
        poly = Polynomial.from_expression('2 - x0 - x2*x0**8', ['x0','x1','x2'])
        singular_set = find_singular_sets_at_one(poly)
        self.assertEqual(singular_set, [(0,2)])

#@attr('active')
class TestSplit(unittest.TestCase):
    seed = 138

    #@attr('active')
    def test_nontrivial_Jacobian(self):
        sector = Sector(
                           cast     = [Polynomial.from_expression('1 -   x0 + x1', ['x0','x1'])],
                           Jacobian =  Polynomial.from_expression('3 - 4*x0     ', ['x0','x1'])
                       )

        split_sectors = list(  split(sector, self.seed+1, 0)  )

        # split at 1/4

        target_split_Jacobian_0 = sympify_expression('(3 - 4 *      x0 *    1/4  ) *    1/4 ')
        target_split_Jacobian_1 = sympify_expression('(3 - 4 * (1 - x0 * (1-1/4))) * (1-1/4)')

        target_split_cast_0 = sympify_expression('1 -      x0 *    1/4   + x1')
        target_split_cast_1 = sympify_expression('1 - (1 - x0 * (1-1/4)) + x1')

        self.assertEqual(len(split_sectors), 2)

        self.assertEqual( (  sympify_expression(split_sectors[0].cast[0]) - target_split_cast_0  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(split_sectors[0].Jacobian) - target_split_Jacobian_0 ).simplify() , 0)

        self.assertEqual( (  sympify_expression(split_sectors[1].cast[0]) - target_split_cast_1  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(split_sectors[1].Jacobian) - target_split_Jacobian_1 ).simplify() , 0)

    #@attr('active')
    def test_splitting(self):
        polysymbols = ['x0','x1','x2']

        cast_poly = Polynomial.from_expression('1 - x0 + x1 - x2', polysymbols)
        other_poly_1 = Polynomial.from_expression('x0**2 * x2**5', polysymbols)
        other_poly_2 = ExponentiatedPolynomial([[0,1,0],[0,0,2]], ['a', 'b'], sp.symbols('other_2_exponent'), polysymbols)
        sector = Sector([cast_poly], [other_poly_1, other_poly_2])

        subsectors = list(  split(sector, self.seed, 0, 2)  )

        # split at [1/2,9/10]

        target_split_Jacobian_0 = sympify_expression('   1/2  *    9/10 ')
        target_split_Jacobian_1 = sympify_expression('   1/2  * (1-9/10)')
        target_split_Jacobian_2 = sympify_expression('(1-1/2) *    9/10 ')
        target_split_Jacobian_3 = sympify_expression('(1-1/2) * (1-9/10)')

        target_split_cast_0 = sympify_expression('1 - x0/2 + x1 - x2*9/10')
        target_split_cast_1 = sympify_expression('1 - x0/2 + x1 - (1-x2*(1-9/10))')
        target_split_cast_2 = sympify_expression('1 - (1-x0/2) + x1 - x2*9/10')
        target_split_cast_3 = sympify_expression('1 - (1-x0/2) + x1 - (1-x2*(1-9/10))')

        target_split_other_1_0 = sympify_expression('(  x0/2)**2 * (  x2*   9/10 )**5')
        target_split_other_1_1 = sympify_expression('(  x0/2)**2 * (1-x2*(1-9/10))**5')
        target_split_other_1_2 = sympify_expression('(1-x0/2)**2 * (  x2*   9/10 )**5')
        target_split_other_1_3 = sympify_expression('(1-x0/2)**2 * (1-x2*(1-9/10))**5')

        target_split_other_2_0 = sympify_expression('(  a * x1  +  b * (  x2*   9/10) **2  )   **   other_2_exponent')
        target_split_other_2_1 = sympify_expression('(  a * x1  +  b * (1-x2*(1-9/10))**2  )   **   other_2_exponent')
        target_split_other_2_2 = sympify_expression('(  a * x1  +  b * (  x2*   9/10) **2  )   **   other_2_exponent')
        target_split_other_2_3 = sympify_expression('(  a * x1  +  b * (1-x2*(1-9/10))**2  )   **   other_2_exponent')

        self.assertEqual(len(subsectors), 4)

        self.assertEqual( (  sympify_expression(subsectors[0].cast[0]) - target_split_cast_0  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[0].other[0]) - target_split_other_1_0  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[0].other[1]) - target_split_other_2_0  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[0].Jacobian) - target_split_Jacobian_0 ).simplify() , 0)

        self.assertEqual( (  sympify_expression(subsectors[1].cast[0]) - target_split_cast_1  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[1].other[0]) - target_split_other_1_1  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[1].other[1]) - target_split_other_2_1  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[1].Jacobian) - target_split_Jacobian_1 ).simplify() , 0)

        self.assertEqual( (  sympify_expression(subsectors[2].cast[0]) - target_split_cast_2  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[2].other[0]) - target_split_other_1_2  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[2].other[1]) - target_split_other_2_2  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[2].Jacobian) - target_split_Jacobian_2 ).simplify() , 0)

        self.assertEqual( (  sympify_expression(subsectors[3].cast[0]) - target_split_cast_3  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[3].other[0]) - target_split_other_1_3  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[3].other[1]) - target_split_other_2_3  ).simplify() , 0)
        self.assertEqual( (  sympify_expression(subsectors[3].Jacobian) - target_split_Jacobian_3 ).simplify() , 0)

#@attr('active')
class TestSplitSingular(unittest.TestCase):
    seed = 138

    #@attr('active')
    def test_split_singular(self):
        poly = Polynomial.from_expression('1 - x0 + A*x1 - x2', ['x0','x1','x2'])
        initial_sector = Sector([poly])

        split_sectors = list( split_singular(initial_sector, self.seed) )

        # singular at 1 for "x0->1" or "x2->1"; expect split at for "x0" and "x2" (at [1/2,9/10])
        target_splits = []
        target_splits.append(   sympify_expression('1 -    x0/2  + A*x1 -    x2*   9/10  ')   )
        target_splits.append(   sympify_expression('1 -    x0/2  + A*x1 - (1-x2*(1-9/10))')   )
        target_splits.append(   sympify_expression('1 - (1-x0/2) + A*x1 -    x2*   9/10  ')   )
        target_splits.append(   sympify_expression('1 - (1-x0/2) + A*x1 - (1-x2*(1-9/10))')   )

        target_Jacobians = sympify_expression([
                                           '   1/2  *    9/10',
                                           '   1/2  * (1-9/10)',
                                           '(1-1/2) *    9/10 ',
                                           '(1-1/2) * (1-9/10)'
                                      ])

        self.assertEqual(len(split_sectors), 4)
        for i in range(4):
            print(i)
            sympified_split_poly = sympify_expression(split_sectors[i].cast[0])
            self.assertEqual(  (sympified_split_poly - target_splits[i]).simplify() , 0  )

            sympyfied_split_Jacobian = sympify_expression(split_sectors[i].Jacobian)
            self.assertEqual(  (sympyfied_split_Jacobian - target_Jacobians[i]).simplify() , 0  )

    #@attr('active')
    def test_selected_parameters(self):
        polysymbols = ['x0','x1','x2','eps']
        poly = Polynomial.from_expression('1 - x0 + A*x1 - x2', polysymbols)
        exponentiated_poly = (  poly ** Polynomial.from_expression('1 - 2*eps' , polysymbols)  ).simplify()
        initial_sector = Sector([exponentiated_poly])

        split_sectors = list( split_singular(initial_sector, self.seed, indices=[0,1,2]) )

        # singular at 1 for "x0->1" or "x2->1"; expect split at for "x0" and "x2" (at [1/2,9/10])
        target_splits = []
        target_splits.append(   sympify_expression('(1 -    x0/2  + A*x1 -    x2*   9/10  ) ** (1 - 2*eps)')   )
        target_splits.append(   sympify_expression('(1 -    x0/2  + A*x1 - (1-x2*(1-9/10))) ** (1 - 2*eps)')   )
        target_splits.append(   sympify_expression('(1 - (1-x0/2) + A*x1 -    x2*   9/10  ) ** (1 - 2*eps)')   )
        target_splits.append(   sympify_expression('(1 - (1-x0/2) + A*x1 - (1-x2*(1-9/10))) ** (1 - 2*eps)')   )

        target_Jacobian = sympify_expression('1/4')

        self.assertEqual(len(split_sectors), 4)
        for i in range(4):
            print(i)
            sympified_split_poly = sympify_expression(split_sectors[i].cast[0])
            self.assertEqual(  (sympified_split_poly - target_splits[i]).simplify() , 0  )
