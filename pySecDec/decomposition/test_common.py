"""Unit tests for the Sector container class"""

from .common import *
from .common import _sector2array, _collision_safe_hash
from ..algebra import Polynomial, ExponentiatedPolynomial, Product
from ..matrix_sort import iterative_sort, Pak_sort, light_Pak_sort
import unittest
import sympy as sp
import numpy as np
from itertools import permutations
from nose.plugins.attrib import attr
import sys

python_major_version = sys.version[0]

try:
    # use "$SECDEC_CONTRIB/bin/dreadnaut" if "$SECDEC_CONTRIB" is defined
    dreadnaut_executable = os.path.join(os.environ['SECDEC_CONTRIB'], 'bin', 'dreadnaut')
except KeyError:
    # "$SECDEC_CONTRIB" is not defined --> let the system find "dreadnaut"
    dreadnaut_executable = 'dreadnaut'

class TestSector(unittest.TestCase):
    def setUp(self):
        self.poly = Polynomial([(1,0,0,4),(0,1,0,1),(0,0,1,0)],[1,1,1])
        self.sector = Sector([self.poly])

    #@attr('active')
    def test_init(self):
        # Feynman parameters are the ti
        # input is part of the 1Loop box

        # F = -s12*t1 - s23*t0*t2
        F = Polynomial([(0,1,0),(1,0,1)],["-s12","-s23"])

        # U = 1 + t0 + t1 + t2
        U = Polynomial([(0,0,0),(1,0,0),(0,1,0),(0,0,1)],[1,1,1,1])

        # "empty" Jacobian in the sense that it is
        # the constant Polynomial with unit constant
        Jacobian = Polynomial([(0,0,0)],[1])

        other_polynomial = Polynomial([(0,1,2),(1,0,5),(1,2,3),(9,4,2)],[1,'A','C','g'], polysymbols=['x','y','z'])

        self.assertRaisesRegexp(AssertionError, 'number of variables.*equal', Sector, [F], [self.poly])
        self.assertRaisesRegexp(AssertionError, '(f|F)irst factor.*monomial', Sector, [Product(F,U)])
        self.assertRaisesRegexp(AssertionError, 'two factors', Sector, [Product(F,U,Jacobian)])
        self.assertRaisesRegexp(AssertionError, 'at least one', Sector, [])
        self.assertRaisesRegexp(AssertionError, 'other.*type.*Polynomial', Sector, [F], other=[Product(Jacobian,U)])
        Sector([Product(Jacobian,F)])

        sector = Sector([F])
        self.assertEqual(str(sector.Jacobian), str(Jacobian))

        sector = Sector([other_polynomial]) # constructor should factorize
        self.assertEqual(str(sector.cast[0]), '( + (1)*z**2) * ( + (1)*y + (A)*x*z**3 + (C)*x*y**2*z + (g)*x**9*y**4)')

    def test_keep_exponent(self):
        exponentiated_poly = ExponentiatedPolynomial(self.poly.expolist, self.poly.coeffs, polysymbols=self.poly.polysymbols, exponent='4-2*eps')
        sector = Sector([exponentiated_poly])
        for i in range(2):
            self.assertTrue(type(sector.cast[0].factors[i]) is ExponentiatedPolynomial)
            self.assertEqual(  (sector.cast[0].factors[i].exponent - sp.sympify('4-2*eps')).simplify() , 0  )

    def test_access(self):
        self.assertEqual(self.sector.other,[])
        self.assertEqual(len(self.sector.cast),1)
        self.assertEqual(str(self.sector.cast[0].factors[1]),str(self.poly))

    def test_copy(self):
        sector = self.sector.copy()
        self.assertEqual(sector.other,self.sector.other)
        self.assertEqual(len(self.sector.cast),len(sector.cast))
        self.assertEqual(str(self.sector.cast[0].factors[1]),str(sector.cast[0].factors[1]))
        self.assertEqual(self.sector.number_of_variables,sector.number_of_variables)

        # really made a copy?
        sector.cast[0].factors[1].expolist += 1
        self.assertNotEqual(str(self.sector.cast[0].factors[1]),sector.cast[0].factors[1])

#@attr('active')
class TestSymmetryFinding(unittest.TestCase):
    def setUp(self):
        self.Jacobian = Polynomial([(1,0)], ['a'])

        self.p0 = Polynomial([(0,1),(1,3)], ['a','b'])

        self.p1_mono = Polynomial([(1,0)], [1])
        self.p1_poly = Polynomial([(2,2),(2,1)], ['c','d'])
        self.p1 = Product(self.p1_mono, self.p1_poly)

        self.p2 = Polynomial([(2,0),(0,1),(1,1)], ['e','f','g'])

        self.sector = Sector([self.p0,self.p1], [self.p2], self.Jacobian)


        self.swapped_p0 = Polynomial([(1,0),(3,1)], ['a','b'])
        self.swapped_Jacobian = Polynomial([(0,1)], ['swapped_Jacobian_coeff'])
        self.sector_p0 = Sector([self.p0], Jacobian=self.Jacobian)
        self.sector_swapped_p0 = Sector([self.swapped_p0], Jacobian=self.swapped_Jacobian)

        # hard example: from "examples/triangle"
        self.symbols_hard = ['x%i'%i for i in range(6)]
        self.Jacobian_hard = Polynomial([[1]*len(self.symbols_hard)], ['a'])
        self.hard_p1 = (Polynomial.from_expression(''' + (1)*x1*x2 + (1)*x1*x3 + (1)*x1*x4 + (1)*x1*x5 + (1)*x1 + (1)*x2*x3 + (1)*x2*x4
                                                       + (1)*x2*x5 + (1)*x3 + (1)*x4 + (1)*x5 ''', self.symbols_hard) ** sp.sympify('3*eps')
                       ).simplify()
        self.hard_p2 = (Polynomial.from_expression(''' + (s)*x1**2*x2 + (s)*x1**2*x3 + (s)*x1**2*x4 + (s)*x1**2*x5 + (s)*x1**2
                                                       + (s)*x1*x2**2 + (2*s)*x1*x2*x3 + (s)*x1*x2*x4 + (2*s)*x1*x2*x5 + (s)*x1*x2
                                                       + (2*s)*x1*x3 + (-s)*x1*x4*x5 + (2*s)*x1*x4 + (s)*x1*x5 + (s)*x1
                                                       + (s)*x2**2*x3 + (s)*x2**2*x4 + (s)*x2**2*x5 + (s)*x2*x3
                                                       + (-s)*x2*x4*x5 + (s)*x2*x4 + (s)*x2*x5 + (s)*x3 + (-s)*x4*x5
                                                       + (s)*x4 + (s)*x5 ''', self.symbols_hard) ** sp.sympify('-2*eps - 2')
                       ).simplify()
        self.hard_p1_permuted = (Polynomial.from_expression(''' + (1)*x0*x1 + (1)*x0*x3 + (1)*x0*x4
                                                                + (1)*x0*x5 + (1)*x1*x3 + (1)*x1*x4 + (1)*x1*x5 + (1)*x1
                                                                + (1)*x3 + (1)*x4 + (1)*x5''', self.symbols_hard) ** sp.sympify('3*eps')
                                ).simplify()
        self.hard_p2_permuted = (Polynomial.from_expression(''' + (s)*x0**2*x1 + (s)*x0**2*x3 + (s)*x0**2*x4 + (s)*x0**2*x5 + (s)*x0*x1**2
                                                                + (2*s)*x0*x1*x3 + (2*s)*x0*x1*x4 + (s)*x0*x1*x5 + (s)*x0*x1 + (s)*x0*x3
                                                                + (-s)*x0*x4*x5 + (s)*x0*x4 + (s)*x0*x5 + (s)*x1**2*x3 + (s)*x1**2*x4
                                                                + (s)*x1**2*x5 + (s)*x1**2 + (2*s)*x1*x3 + (-s)*x1*x4*x5 + (s)*x1*x4
                                                                + (2*s)*x1*x5 + (s)*x1 + (s)*x3 + (-s)*x4*x5 + (s)*x4
                                                                + (s)*x5''', self.symbols_hard) ** sp.sympify('-2*eps - 2')
                                ).simplify()

        self.sector_hard = Sector([self.hard_p1,self.hard_p2], [], self.Jacobian_hard)
        self.sector_swapped_hard = Sector([self.hard_p1_permuted,self.hard_p2_permuted], [], self.Jacobian_hard)

        self.a, self.b, self.c, self.d, self.e, self.f, self.g = sp.symbols('a b c d e f g')

    #@attr('active')
    def test_sector2array(self):
        SecDecInternalCast, SecDecInternalOther = sp.symbols('SecDecInternalCast SecDecInternalOther')

        combined_expolists, combined_coeffs = _sector2array(self.sector)

        # note that `Sector` factorizes on construction
        target_combined_expolists = np.array([
                                                (1,0),            # Jacobian
                                                (0,1),(1,3),      # p0
                                                (3,1),(3,2),      # p1 --> note the reordering
                                                (2,0),(0,1),(1,1) # p2
                                            ])
        target_combined_coeffs = np.array([
                                               # Jacobian coefficient is ignored (replaced by a dummy)
                                               1,

                                               # p0
                                               self.a*SecDecInternalCast(0),self.b*SecDecInternalCast(0),

                                               # p1
                                               self.d*SecDecInternalCast(1),self.c*SecDecInternalCast(1),

                                               # p2
                                               self.e*SecDecInternalOther(0),self.f*SecDecInternalOther(0),self.g*SecDecInternalOther(0)
                                         ])

        np.testing.assert_array_equal(combined_expolists, target_combined_expolists)
        np.testing.assert_array_equal(combined_coeffs, target_combined_coeffs)

    #@attr('active')
    def test_collision_safe_hash(self):
        class CustomHash(object):
            def __init__(self, hash, value):
                self.hash = hash
                self.value = value
            def __hash__(self):
                return self.hash
            def __eq__(self, other):
                if isinstance(other, CustomHash):
                    return self.value == other.value
                else:
                    return NotImplemented
            def __ne__(self, other):
                if isinstance(other, CustomHash):
                    return self.value != other.value
                else:
                    return NotImplemented
            def __str__(self):
                return "CustomHash(hash=%i,value=%i)" % (self.hash,self.value)
            __repr__ = __str__

        array_with_hash_collisions = np.array([
            CustomHash(1,1), CustomHash(2,2), CustomHash(1,1), CustomHash(1,4), CustomHash(2,5)
        ]) # hash collision since ``CustomHash(n,1) != CustomHash(n,2)`` but hashes are equal

        array_collisions_resolved = _collision_safe_hash(array_with_hash_collisions)

        for i,j in permutations(range(len(array_collisions_resolved)), 2):
            print(i,j)
            if i==0 and j==2 or i==2 and j==0:
                self.assertEqual(array_collisions_resolved[i], array_collisions_resolved[j])
            else:
                self.assertNotEqual(array_collisions_resolved[i], array_collisions_resolved[j])

    #@attr('active')
    def test_squash_symmetry_redundant_sectors_2D(self):
        sectors = [self.sector_p0.copy(), self.sector_swapped_p0.copy()]

        # test symmetry finding by sorting
        for sort_function in (iterative_sort, Pak_sort):
            reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors, sort_function)

            self.assertEqual(len(reduced_sectors), 1)
            self.assertEqual(reduced_sectors[0].Jacobian.coeffs[0], sp.sympify('a+swapped_Jacobian_coeff'))
            self.assertEqual( (sp.sympify(reduced_sectors[0].cast[0]) - sp.sympify(self.p0.copy())).simplify() , 0 )

        # test symmetry finding by graph (using dreadnaut)
        reduced_sectors = squash_symmetry_redundant_sectors_dreadnaut(sectors, dreadnaut_executable, workdir='tmpdir_test_squash_symmetry_redundant_sectors_2D_python' + python_major_version)
        self.assertEqual(len(reduced_sectors), 1)
        self.assertEqual(reduced_sectors[0].Jacobian.coeffs[0], sp.sympify('a+swapped_Jacobian_coeff'))
        self.assertEqual((sp.sympify(reduced_sectors[0].cast[0]) - sp.sympify(self.p0.copy())).simplify(), 0)

    #@attr('active')
    def test_squash_symmetry_hard(self):
        sectors = [self.sector_hard.copy(), self.sector_swapped_hard.copy()]

        # test symmetry finding by iterative sorting, fails
        reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors, iterative_sort)
        self.assertNotEqual(len(reduced_sectors), 1)
        self.assertNotEqual(reduced_sectors[0].Jacobian.coeffs[0], sp.sympify('2*a'))

        # test symmetry finding using light sorting, fails
        reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors, light_Pak_sort)
        self.assertNotEqual(len(reduced_sectors), 1)
        self.assertNotEqual(reduced_sectors[0].Jacobian.coeffs[0], sp.sympify('2*a'))

        # test symmetry finding by graph (using dreadnaut)
        reduced_sectors = squash_symmetry_redundant_sectors_dreadnaut(sectors, dreadnaut_executable, workdir='tmpdir_test_squash_symmetry_hard_python' + python_major_version)
        self.assertEqual(len(reduced_sectors), 1)
        self.assertEqual(reduced_sectors[0].Jacobian.coeffs[0], sp.sympify('2*a'))

        # test symmetry finding using Pak's full algorithm
        reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors, Pak_sort)
        self.assertEqual(len(reduced_sectors), 1)
        self.assertEqual(reduced_sectors[0].Jacobian.coeffs[0], sp.sympify('2*a'))

    #@attr('active')
    def test_symmetry_4D(self):
        # sectors 0 and 2 are related by permutation, sector 1 is unrelated
        sector0_p0 = Polynomial([(0,1,1,3),(2,2,4,3)], ['a','b'])
        sector0_p1 = Polynomial([(1,2,1,2),(1,2,3,1)], ['1','1'])
        sector0 = Sector([sector0_p0, sector0_p1])

        sector1_p0 = Polynomial([(0,5,1,3),(2,2,4,3)], ['a','b'])
        sector1_p1 = Polynomial([(1,2,2,1),(3,2,1,1)], ['1','1'])
        sector1 = Sector([sector1_p0, sector1_p1])

        sector2_p0 = Polynomial([(4,2,3,2),(1,1,3,0)], ['b','a'])
        sector2_p1 = Polynomial([(1,2,2,1),(3,2,1,1)], ['1','1'])
        sector2 = Sector([sector2_p0, sector2_p1])

        for i in range(2): # run twice to check if the variables `sectorI` are not modified
            sectors_with_redundancy = (sector0, sector1, sector2)

            # test symmetry finding by sorting
            for sort_function in (iterative_sort, Pak_sort):
                reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors_with_redundancy, sort_function)

                # should have found the symmetry and pruned `sector0` or `sector2`
                self.assertEqual(len(reduced_sectors), 2)

                # should have either `sector0` or `sector2` in `reduced_sectors`
                have_sector_0 = (str(reduced_sectors[0].cast) == str(sector0.cast) or str(reduced_sectors[1].cast) == str(sector0.cast))
                self.assertTrue( (str(reduced_sectors[0].cast) == str(sector0.cast if have_sector_0 else sector2.cast))
                              or (str(reduced_sectors[1].cast) == str(sector0.cast if have_sector_0 else sector2.cast)) )

                # `sector1` should be untouched and Jacobian coefficient should have been increased by one
                self.assertTrue( (str(reduced_sectors[0].Jacobian) == ' + (2)' and str(reduced_sectors[1]) == str(sector1))
                              or (str(reduced_sectors[1].Jacobian) == ' + (2)' and str(reduced_sectors[0]) == str(sector1)) )

            # test symmetry finding by graph (using dreadnaut)
            reduced_sectors=squash_symmetry_redundant_sectors_dreadnaut(sectors_with_redundancy, dreadnaut_executable, workdir='tmpdir_test_symmetry_4D_python' + python_major_version)

            # should have found the symmetry and pruned `sector0` or `sector2`
            self.assertEqual(len(reduced_sectors), 2)

            # should have either `sector0` or `sector2` in `reduced_sectors`
            have_sector_0 = (str(reduced_sectors[0].cast) == str(sector0.cast) or str(reduced_sectors[1].cast) == str(sector0.cast))
            self.assertTrue((str(reduced_sectors[0].cast) == str(sector0.cast if have_sector_0 else sector2.cast))
                        or (str(reduced_sectors[1].cast) == str(sector0.cast if have_sector_0 else sector2.cast)))

            # `sector1` should be untouched and Jacobian coefficient should have been increased by one
            self.assertTrue((str(reduced_sectors[0].Jacobian) == ' + (2)' and str(reduced_sectors[1]) == str(sector1))
                        or (str(reduced_sectors[1].Jacobian) == ' + (2)' and str(reduced_sectors[0]) == str(sector1)))

    #@attr('active')
    def test_symmetry_special_sorting(self):
        # sectors 0 and 1 are related by permutation
        sector0_p0 = Polynomial([(0,1,2,3),(3,2,1,0)], ['a','a'])
        sector0_p1 = Polynomial([(1,2,1,2),(1,2,3,1)], ['1','1'])
        sector0 = Sector([sector0_p0, sector0_p1])

        sector1_p0 = Polynomial([(3,1,2,0),(0,2,1,3)], ['a','a'])
        sector1_p1 = Polynomial([(1,3,2,1),(1,1,2,2)], ['1','1'])
        sector1 = Sector([sector1_p0, sector1_p1])

        for i in range(2): # run twice to check if the variables `sectorI` are not modified
            sectors_with_redundancy = (sector0, sector1)

            for sort_function in (iterative_sort, Pak_sort):
                reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors_with_redundancy, sort_function)

                # should have found the symmetry
                self.assertEqual(len(reduced_sectors), 1)

                # should have either `sector0` or `sector1` in `reduced_sectors` with Jacobian doubled
                have_sector_0 = (str(reduced_sectors[0].cast) == str(sector0.cast) or str(reduced_sectors[1].cast) == str(sector0.cast))
                target_reduced_sectors = [sector0.copy() if have_sector_0 else sector1.copy()]
                target_reduced_sectors[0].Jacobian.coeffs[0] = 2
                self.assertEqual( str(reduced_sectors), str(target_reduced_sectors) )

            reduced_sectors = squash_symmetry_redundant_sectors_dreadnaut(sectors_with_redundancy, dreadnaut_executable, workdir='tmpdir_test_symmetry_special_sorting_python' + python_major_version)

            # should have found the symmetry
            self.assertEqual(len(reduced_sectors), 1)

            # should have either `sector0` or `sector1` in `reduced_sectors` with Jacobian doubled
            have_sector_0 = (str(reduced_sectors[0].cast) == str(sector0.cast) or str(reduced_sectors[1].cast) == str(sector0.cast))
            target_reduced_sectors = [sector0.copy() if have_sector_0 else sector1.copy()]
            target_reduced_sectors[0].Jacobian.coeffs[0] = 2
            self.assertEqual(str(reduced_sectors), str(target_reduced_sectors))

    #@attr('active')
    def test_symmetry_same_term_in_different_polynomials(self):
        eps = sp.symbols('eps')

        # sectors 0 and 2 are related by permutation, sector 1 is unrelated (due to exponent)
        sector0_p0 = Polynomial([(0,1,1,3),(2,2,4,3)], ['a','b'])
        sector0_p1 = Polynomial([(1,2,1,2),(1,2,3,1)], ['1','1'])
        sector0 = Sector([sector0_p0, sector0_p1])

        sector1_p0 = ExponentiatedPolynomial([(0,1,1,3),(2,2,4,3)], ['a','b'], exponent=eps)
        sector1_p1 = ExponentiatedPolynomial([(1,2,1,2),(1,2,3,1)], ['1','1'], exponent=eps)
        sector1 = Sector([sector1_p0, sector1_p1])

        sector2_p0 = Polynomial([(4,2,3,2),(1,1,3,0)], ['b','a'])
        sector2_p1 = Polynomial([(1,2,2,1),(3,2,1,1)], ['1','1'])
        sector2 = Sector([sector2_p0, sector2_p1])

        for i in range(2): # run twice to check if the variables `sectorI` are not modified
            print(i)
            sectors_with_redundancy = (sector0, sector1, sector2)

            for sort_function in (iterative_sort, Pak_sort):
                reduced_sectors = squash_symmetry_redundant_sectors_sort(sectors_with_redundancy, sort_function)

                # should have found the symmetry and pruned `sector0` or `sector2`
                self.assertEqual(len(reduced_sectors), 2)

                # should have either `sector0` or `sector2` in `reduced_sectors`
                have_sector_0 = (str(reduced_sectors[0].cast) == str(sector0.cast) or str(reduced_sectors[1].cast) == str(sector0.cast))
                self.assertTrue(str(reduced_sectors[0].cast) == str(sector0.cast if have_sector_0 else sector2.cast)
                             or str(reduced_sectors[1].cast) == str(sector0.cast if have_sector_0 else sector2.cast))

                # Jacobian coefficient should have been increased by one while `sector1` should be untouched
                self.assertTrue( (str(reduced_sectors[0].Jacobian) == ' + (2)' and str(reduced_sectors[1]) == str(sector1))
                              or (str(reduced_sectors[1].Jacobian) == ' + (2)' and str(reduced_sectors[0]) == str(sector1)) )

            reduced_sectors = squash_symmetry_redundant_sectors_dreadnaut(sectors_with_redundancy, dreadnaut_executable, workdir='tmpdir_test_symmetry_same_term_in_different_polynomials_python' + python_major_version)

            # should have found the symmetry and pruned `sector0` or `sector2`
            self.assertEqual(len(reduced_sectors), 2)

            # should have either `sector0` or `sector2` in `reduced_sectors`
            have_sector_0 = (
            str(reduced_sectors[0].cast) == str(sector0.cast) or str(reduced_sectors[1].cast) == str(sector0.cast))
            self.assertTrue(str(reduced_sectors[0].cast) == str(sector0.cast if have_sector_0 else sector2.cast)
                            or str(reduced_sectors[1].cast) == str(sector0.cast if have_sector_0 else sector2.cast))

            # Jacobian coefficient should have been increased by one while `sector1` should be untouched
            self.assertTrue((str(reduced_sectors[0].Jacobian) == ' + (2)' and str(reduced_sectors[1]) == str(sector1))
                            or (
                            str(reduced_sectors[1].Jacobian) == ' + (2)' and str(reduced_sectors[0]) == str(sector1)))


class TestOther(unittest.TestCase):
    def test_refactorize(self):
        prod = Product(Polynomial([(0,0,0)],[1],'t'), Polynomial([(1,1,0),(1,0,1)],["-s12","-s23"],'t'))

        self.assertEqual(str(prod.factors[0]), ' + (1)')
        self.assertEqual(str(prod.factors[1]), ' + (-s12)*t0*t1 + (-s23)*t0*t2')

        copy0 = prod.copy()
        copy1 = prod.copy()
        copy2 = prod.copy()

        refactorize(copy0) # refactorize all parameters -> should find the factorization of parameter 0
        refactorize(copy1,0) # refactorize parameter 0 -> should find a factorization
        refactorize(copy2,1) # refactorize parameter 1 -> should NOT find the factorization of parameter 0

        self.assertEqual(str(copy0.factors[0]), str(copy1.factors[0]))
        self.assertEqual(str(copy0.factors[1]), str(copy1.factors[1]))

        self.assertEqual(str(copy1.factors[0]), ' + (1)*t0')
        self.assertEqual(str(copy1.factors[1]), ' + (-s12)*t1 + (-s23)*t2')

        self.assertEqual(str(copy2.factors[0]), ' + (1)')
        self.assertEqual(str(copy2.factors[1]), ' + (-s12)*t0*t1 + (-s23)*t0*t2')
