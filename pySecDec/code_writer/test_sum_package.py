from .sum_package import *
from ..algebra import Polynomial, ExponentiatedPolynomial
from ..misc import sympify_expression
from ..make_package import MakePackage
from ..loop_integral import LoopIntegralFromGraph, LoopPackage
from nose.plugins.attrib import attr
import sys
import numpy as np
import sympy as sp
import unittest

python_major_version = sys.version[0]

#@attr('active')
class TestCoefficient(unittest.TestCase):
    def setUp(self):
        self.regulators = ['eps']
        self.parameters = ['s','t']

        self.numerator = '''
            2*(-1680*s^6 + 96*(-87*s^5 + 16*s^6)*t - 24*(638*s^4 - 242*s^5 + 21*s^6)*t^2 -
            96*(216*s^3 - 218*s^4 + 69*s^5 - 7*s^6)*t^3 +
            24*(-1496*s^2 + 1862*s^3 - 668*s^4 + 78*s^5 - s^6)*t^4 -
            48*(816*s - 856*s^2 + 291*s^3 - 36*s^4 + s^5)*t^5 +
            24*(-640 + 544*s - 168*s^2 + 22*s^3 - s^4)*t^6 -
            (4-2*eps)^2*(2288*s^6 - 8*(-1292*s^5 + 223*s^6)*t + (14384*s^4 - 4816*s^5 + 327*s^6)*
            t^2 + 2*(6656*s^3 - 10012*s^4 + 3819*s^5 - 433*s^6)*t^3 -
            (-33344*s^2 + 51264*s^3 - 20028*s^4 + 2464*s^5 - 35*s^6)*t^4 +
            2*(23168*s - 25744*s^2 + 9088*s^3 - 1165*s^4 + 35*s^5)*t^5 -
            (-19968 + 17216*s - 5424*s^2 + 732*s^3 - 35*s^4)*t^6) +
            2*(4-2*eps)*(1696*s^6 - 208*(-39*s^5 + 7*s^6)*t + (13472*s^4 - 4916*s^5 + 399*s^6)*
            t^2 + 2*(8184*s^3 - 9274*s^4 + 3135*s^5 - 332*s^6)*t^3 -
            (-31520*s^2 + 42224*s^3 - 15650*s^4 + 1866*s^5 - 25*s^6)*t^4 +
            10*(3744*s - 4016*s^2 + 1386*s^3 - 174*s^4 + 5*s^5)*t^5 -
            (-15232 + 13024*s - 4056*s^2 + 538*s^3 - 25*s^4)*t^6) +
            2*(4-2*eps)^3*(320*s^6 - 16*(-85*s^5 + 14*s^6)*t + (1472*s^4 - 412*s^5 + 15*s^6)*t^2 +
            2*(296*s^3 - 1030*s^4 + 471*s^5 - 58*s^6)*t^3 -
            (-3296*s^2 + 6288*s^3 - 2622*s^4 + 334*s^5 - 5*s^6)*t^4 +
            2*(2912*s - 3376*s^2 + 1222*s^3 - 160*s^4 + 5*s^5)*t^5 -
            (-2688 + 2336*s - 744*s^2 + 102*s^3 - 5*s^4)*t^6) +
            (4-2*eps)^4*(-64*s^6 + 8*(-32*s^5 + 5*s^6)*t + (-192*s^4 + 32*s^5 + 3*s^6)*t^2 +
            2*(64*s^3 + 132*s^4 - 81*s^5 + 11*s^6)*t^3 +
            (-384*s^2 + 1072*s^3 - 484*s^4 + 64*s^5 - s^6)*t^4 -
            2*(512*s - 624*s^2 + 232*s^3 - 31*s^4 + s^5)*t^5 +
            (-512 + 448*s - 144*s^2 + 20*s^3 - s^4)*t^6))
        '''
        self.sympified_numerator = sympify_expression(self.numerator)

        self.denominator = '''
            2048*s^6*t - 3072*(-4*s^5 + s^6)*t^2 + 1792*(16*s^4 - 8*s^5 + s^6)*t^3 +
            512*(64*s^3 - 48*s^4 + 12*s^5 - s^6)*t^4 -
            72*(-256*s^2 + 256*s^3 - 96*s^4 + 16*s^5 - s^6)*t^5 +
            4*(1024*s - 1280*s^2 + 640*s^3 - 160*s^4 + 20*s^5 - s^6)*t^6 -
            (4-2*eps)*(512*s^6*t - 768*(-4*s^5 + s^6)*t^2 + 448*(16*s^4 - 8*s^5 + s^6)*t^3 +
            128*(64*s^3 - 48*s^4 + 12*s^5 - s^6)*t^4 -
            18*(-256*s^2 + 256*s^3 - 96*s^4 + 16*s^5 - s^6)*t^5 +
            (1024*s - 1280*s^2 + 640*s^3 - 160*s^4 + 20*s^5 - s^6)*t^6)
        '''
        self.sympified_denominator = sympify_expression(self.denominator)

    #@attr('active')
    @attr('slow')
    def test_coefficient_order_num_den(self):
        coeff = Coefficient([self.numerator],[self.denominator], self.parameters)
        lowest_orders, coefficient_definition = coeff.process(self.regulators, workdir='tmpdir_test_coefficient_order_num_den' + python_major_version)

        # check lowest_orders
        np.testing.assert_array_equal(lowest_orders, [-1])

        # check coefficient_definition
        processed_expressions = coefficient_definition.split(';')
        processed_numerator = sympify_expression( processed_expressions[0].split('=')[1] )
        processed_denominator = sympify_expression( processed_expressions[1].split('=')[1] )
        processed_regulator_factor = sympify_expression( processed_expressions[2].split('=')[1] )
        self.assertEqual(
            (
                (self.sympified_numerator / self.sympified_denominator) / \
                (processed_numerator / processed_denominator * processed_regulator_factor)
            ).simplify()
            , 1
        )

    #@attr('active')
    @attr('slow')
    def test_coefficient_order_empty(self):
        coeff = Coefficient([],[], self.parameters)
        lowest_orders, coefficient_definition = coeff.process(self.regulators, workdir='tmpdir_test_coefficient_order_empty' + python_major_version)

        # check lowest_orders
        np.testing.assert_array_equal(lowest_orders, [0])

        # check coefficient_definition
        processed_expressions = coefficient_definition.split(';')
        processed_numerator = sympify_expression( processed_expressions[0].split('=')[1] )
        processed_denominator = sympify_expression( processed_expressions[1].split('=')[1] )
        processed_regulator_factor = sympify_expression( processed_expressions[2].split('=')[1] )
        for item in (processed_numerator, processed_denominator, processed_regulator_factor):
            self.assertEqual(item, 1)

    #@attr('active')
    @attr('slow')
    def test_coefficient_order_no_numerator(self):
        coeff = Coefficient([],[self.denominator,'eps^2'], self.parameters)
        lowest_orders, coefficient_definition = coeff.process(self.regulators, workdir='tmpdir_test_coefficient_order_no_numerator' + python_major_version)

        # check lowest_orders
        np.testing.assert_array_equal(lowest_orders, [-3])

        # check coefficient_definition
        processed_expressions = coefficient_definition.split(';')
        processed_numerator = sympify_expression( processed_expressions[0].split('=')[1] )
        processed_denominator = sympify_expression( processed_expressions[1].split('=')[1] )
        processed_regulator_factor = sympify_expression( processed_expressions[2].split('=')[1] )
        self.assertEqual(
            (
                (1 / self.sympified_denominator / sympify_expression('eps^2')) / \
                (processed_numerator / processed_denominator * processed_regulator_factor)
            ).simplify()
            , 1
        )

    #@attr('active')
    @attr('slow')
    def test_coefficient_order_no_denominator(self):
        coeff = Coefficient([self.numerator,'eps^2'], [], self.parameters)
        lowest_orders, coefficient_definition = coeff.process(self.regulators, workdir='tmpdir_test_coefficient_order_no_denominator' + python_major_version)

        # check lowest_orders
        np.testing.assert_array_equal(lowest_orders, [2])

        # check coefficient_definition
        processed_expressions = coefficient_definition.split(';')
        processed_numerator = sympify_expression( processed_expressions[0].split('=')[1] )
        processed_denominator = sympify_expression( processed_expressions[1].split('=')[1] )
        processed_regulator_factor = sympify_expression( processed_expressions[2].split('=')[1] )
        self.assertEqual(
            (
                (self.sympified_numerator * sympify_expression('eps^2')) / \
                (processed_numerator / processed_denominator * processed_regulator_factor)
            ).simplify()
            , 1
        )

    #@attr('active')
    def test_coefficient_with_imaginary_unit(self):
        coeff = Coefficient(['I*(5+I*8)'], [], [])
        lowest_orders, coefficient_definition = coeff.process([], workdir='tmpdir_test_coefficient_with_imaginary_unit' + python_major_version)

        # check lowest_orders
        np.testing.assert_array_equal(lowest_orders, [])

        # check coefficient_definition
        processed_expressions = coefficient_definition.split(';')
        processed_numerator = processed_expressions[0].split('=')[1]
        processed_denominator = processed_expressions[1].split('=')[1]
        processed_regulator_factor = processed_expressions[2].split('=')[1]

        self.assertFalse( '^2' in processed_numerator )

        self.assertEqual( sympify_expression(processed_numerator) - sympify_expression('5*I-8') , 0)
        self.assertEqual( sympify_expression(processed_denominator) - sympify_expression('1') , 0)
        self.assertEqual( sympify_expression(processed_regulator_factor) - sympify_expression('1') , 0)

# --------------------------------- mid-level tests ---------------------------------
class TestSumPackage(unittest.TestCase):
    'Base class to define the tearDown method.'
    def tearDown(self):
        for tmpdir in self.tmpdirs:
            try:
                shutil.rmtree(tmpdir)
            except OSError as error:
                if error.errno == 2: # no such file or directory --> this is what we want anyway
                    pass
                else: # reraise error otherwise
                    raise

#@attr('active')
class TestSumPackageCall(TestSumPackage):
    def setUp(self):
        self.tmpdirs = ['tmpdir_test_sum_package_python' + python_major_version + '_0',
                        'tmpdir_test_sum_package_python' + python_major_version + '_1',
                        'tmpdir_test_sum_package_python' + python_major_version + '_2',
                        'tmpdir_test_sum_package_python' + python_major_version + '_3']

        self.make_package_input_one = MakePackage(
                                            name='make_package_1',
                                            integration_variables=['z0','z1','z2'],
                                            ibp_power_goal=-1,
                                            regulators=['eps','alpha'],
                                            requested_orders=[1,2],
                                            polynomials_to_decompose=[1,Polynomial([[0,0,0,0,0,0,0],[1,1,1,0,0,0,0]],['-s','-t'],['z0','z1','z2','eps','alpha','U','F'])],
                                            polynomial_names=['U','F'],
                                            other_polynomials=['U*z1 + F'],
                                            prefactor=1,
                                            remainder_expression='DummyFunction(z0,eps)',
                                            functions=['DummyFunction'],
                                            real_parameters=['s','t'],
                                            complex_parameters=[sympify_expression('msq')],
                                            form_optimization_level=2,
                                            form_work_space='500M',
                                            form_memory_use=None,
                                            form_threads=2,
                                            form_insertion_depth=0,
                                            contour_deformation_polynomial=None,
                                            positive_polynomials=[],
                                            decomposition_method='iterative_no_primary'
                                        )
        self.make_package_input_two = MakePackage(
                                            name='make_package_2',
                                            integration_variables=['z0','z1'],
                                            ibp_power_goal=-1,
                                            regulators=['eps','alpha'],
                                            requested_orders=[1,2],
                                            polynomials_to_decompose=[1,Polynomial([[0,0,0,0,0,0],[1,1,0,0,0,0]],['-s','-t'],['z0','z1','eps','alpha','U','F'])],
                                            polynomial_names=['U','F'],
                                            other_polynomials=['U*z0 + F'],
                                            prefactor=1,
                                            remainder_expression='DummyFunction(z0,eps)',
                                            functions=['DummyFunction'],
                                            real_parameters=['s','t'],
                                            complex_parameters=[sympify_expression('msq')],
                                            form_optimization_level=2,
                                            form_work_space='500M',
                                            form_memory_use=None,
                                            form_threads=2,
                                            form_insertion_depth=0,
                                            contour_deformation_polynomial=None,
                                            positive_polynomials=[],
                                            decomposition_method='iterative_no_primary'
                                        )
        self.make_package_input_three = MakePackage(
                                            name='make_package_3',
                                            integration_variables=['z0','z1'],
                                            regulators=['eps'],
                                            requested_orders=[0],
                                            polynomials_to_decompose=[1,Polynomial([[0,0,0],[1,1,0]],['-s','-t'],['z0','z1','eps'])],
                                            real_parameters=['s','t'],
                                            form_optimization_level=2,
                                            form_work_space='500M',
                                            form_memory_use=None,
                                            form_threads=2,
                                            form_insertion_depth=0,
                                            decomposition_method='iterative_no_primary'
                                        )

        li1 = LoopIntegralFromGraph(
                  internal_lines=[[0, [1, 2]], [0, [2, 1]]],
                  external_lines=[['p1', 1], ['p2', 2]],
                  replacement_rules=[('p1*p1', 's'), ('p2*p2', 's'), ('p1*p2', 's')]
              )

        self.loop_package_input_one = LoopPackage(
                                          name='loop_package_1',
                                          loop_integral=li1,
                                          requested_orders=[0]
                                      )

    #@attr('active')
    def test_sum_package_one_make_package(self):
        result = sum_package(
            self.tmpdirs[0],
            [self.make_package_input_one],
            regulators=['eps','alpha'],
            requested_orders=[1,2],
            real_parameters=['s','t'],
            complex_parameters=[sympify_expression('msq')],
            coefficients=None,
            pylink_qmc_transforms=['korobov3x3']
        )
        self.assertEqual(result['name'], 'make_package_1') # Should match name of last package_generator

    #@attr('active')
    def test_sum_package_two_make_package(self):
        result = sum_package(
            self.tmpdirs[1],
            [self.make_package_input_one, self.make_package_input_two],
            regulators=['eps','alpha'],
            requested_orders=[1,2],
            real_parameters=['s','t'],
            complex_parameters=[sympify_expression('msq')],
            coefficients=None,
            pylink_qmc_transforms=['korobov3x3']
        )
        self.assertEqual(result['name'], 'make_package_2') # Should match name of last package_generator

    #@attr('active')
    def test_sum_package_one_loop_package(self):
        result = sum_package(
            self.tmpdirs[2],
            [self.loop_package_input_one],
            regulators=['eps'],
            requested_orders=[0],
            real_parameters=['s'],
            coefficients=None,
            pylink_qmc_transforms=['korobov3x3']
        )
        self.assertEqual(result['name'], 'loop_package_1')  # Should match name of last package_generator

    #@attr('active')
    def test_sum_package_mixed_package(self):
        result = sum_package(
            self.tmpdirs[3],
            [self.loop_package_input_one, self.make_package_input_three],
            regulators=['eps'],
            requested_orders=[0],
            real_parameters=['s','t'],
            coefficients=None,
            pylink_qmc_transforms=['korobov3x3']
        )
        self.assertEqual(result['name'], 'make_package_3')  # Should match name of last package_generator
