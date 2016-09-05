"""Unit tests for the template parser module"""

from __future__ import print_function
from .make_package import *
from .make_package import _convert_input, _make_FORM_definition, \
                          _make_FORM_function_definition, _make_FORM_list, \
                          _derivative_muliindex_to_name, _make_FORM_shifted_orders, \
                          _make_FORM_Series_initilization, _validate, \
                          _make_prefactor_function, _make_CXX_function_declaration
from ..algebra import Function, Polynomial, Product, ProductRule, Sum
from nose.plugins.attrib import attr
import sys, shutil
import unittest

python_major_version = sys.version[0]

class TestMakePackage(unittest.TestCase):
    'Base class to define the tearDown method.'
    def tearDown(self):
        try:
            shutil.rmtree(self.tmpdir)
        except OSError as error:
            if error.errno == 2: # no such file or directory --> this is what we want anyway
                pass
            else: # reraise error otherwise
                raise

# ----------------------------------- parse input -----------------------------------
class TestConvertInput(TestMakePackage):
    def setUp(self):
        self.tmpdir = 'tmpdir_test_convert_input_python' + python_major_version
        self.correct_input = dict(
                                      name=self.tmpdir,
                                      integration_variables=['z0','z1','z2'],
                                      regulators=['eps','alpha'],
                                      requested_orders=[1,2],
                                      polynomials_to_decompose=[1,Polynomial([[0,0,0],[1,1,1]],['-s','-t'],['z0','z1','z2','eps','alpha'])],
                                      polynomial_names=['U','F'],
                                      other_polynomials=['U*z1 + F'],
                                      prefactor=1,
                                      remainder_expression='DummyFunction(z0,eps)',
                                      functions=['DummyFunction'],
                                      real_parameters=['s','t'],
                                      complex_parameters=[sp.sympify('msq')],
                                      form_optimization_level=2,
                                      form_work_space='500M',
                                      form_insertion_depth=0,
                                      stabilize=False,
                                      contour_deformation_polynomial=None,
                                      decomposition_method='iterative_no_primary'
                                 )

    #@attr('active')
    def test_convert_input(self):
        _convert_input(**self.correct_input) # should be ok

        requested_orders_wrong_shape = self.correct_input.copy()
        requested_orders_wrong_shape['requested_orders'] = [[1,1],[0,0]]
        self.assertRaisesRegexp(AssertionError, r'requested_orders.*wrong shape.*is \(2, 2\).*should be \(2,\)', _convert_input, **requested_orders_wrong_shape)

        requested_orders_wrong_length = self.correct_input.copy()
        requested_orders_wrong_length['requested_orders'] = [1,1,1]
        self.assertRaisesRegexp(AssertionError, 'length.*requested_orders.*match.*length.*regulators', _convert_input, **requested_orders_wrong_length)

        polynomials_to_decompose_unrelated_polysymbols = self.correct_input.copy()
        polynomials_to_decompose_unrelated_polysymbols['polynomials_to_decompose'] = ['1', Polynomial([[0,0,0],[1,1,1]],['-s','-t'],['x0','x1','x2'])]
        self.assertRaisesRegexp(ValueError, r"\(\-s\) \+ \(\-t\)\*x0\*x1\*x2.*polynomials_to_decompose.*symbols.*\(is.*x0, x1, x2.*should.*\z0, z1, z2, eps, alpha", _convert_input, **polynomials_to_decompose_unrelated_polysymbols)

        polynomials_to_decompose_wrong_polysymbols_in_exponent = self.correct_input.copy()
        polynomials_to_decompose_wrong_polysymbols_in_exponent['polynomials_to_decompose'] = ['1', ExponentiatedPolynomial([[0,0,0,0,0],[1,1,1,0,0]],['-s','-t'],polysymbols=['z0','z1','z2','eps','alpha'],exponent=Polynomial([[0,0,1]],[1],['x0','x1','x2']))]
        _convert_input(**polynomials_to_decompose_wrong_polysymbols_in_exponent) # should be ok

        polynomials_to_decompose_wrong_polysymbols_in_coeff = self.correct_input.copy()
        polynomials_to_decompose_wrong_polysymbols_in_coeff['polynomials_to_decompose'] = [1,Polynomial([[0,0,0],[1,1,1]],[Polynomial([[0,0,0]], ['-s'], ['eps','alpha']),'-t'],['z0','z1','z2','eps','alpha'])]
        _convert_input(**polynomials_to_decompose_wrong_polysymbols_in_coeff) # should be ok

        polynomials_to_decompose_polynomial_in_coeff = self.correct_input.copy()
        polynomials_to_decompose_polynomial_in_coeff['polynomials_to_decompose'] = [1,Polynomial([[0,0,0],[1,1,1]],[Polynomial([[0,0,0]], ['-s'], ['z0','z1','z2','eps','alpha']),'-t'],['z0','z1','z2','eps','alpha'])]
        _convert_input(**polynomials_to_decompose_polynomial_in_coeff) # should be ok

        polynomials_to_decompose_sympy_exponent = self.correct_input.copy()
        polynomials_to_decompose_sympy_exponent['polynomials_to_decompose'] = ['1', ExponentiatedPolynomial([[0,0,0,0,0],[1,1,1,0,0]],['-s','-t'],polysymbols=['z0','z1','z2','eps','alpha'],exponent='1+eps')]
        _convert_input(**polynomials_to_decompose_sympy_exponent) # should be ok

        polynomials_to_decompose_nontrivial_as_string = self.correct_input.copy()
        polynomials_to_decompose_nontrivial_as_string['polynomials_to_decompose'] = ['1', '(-s -t*z0*z1*z2)**(2-4*eps+alpha)']
        _convert_input(**polynomials_to_decompose_nontrivial_as_string) # should be ok

        polynomials_to_decompose_negative_insertion_depth = self.correct_input.copy()
        polynomials_to_decompose_negative_insertion_depth['form_insertion_depth'] = -3
        self.assertRaisesRegexp(AssertionError, 'form_insertion_depth.*negative', _convert_input, **polynomials_to_decompose_negative_insertion_depth)

        polynomials_to_decompose_noninteger_insertion_depth = self.correct_input.copy()
        polynomials_to_decompose_noninteger_insertion_depth['form_insertion_depth'] = 1.2
        self.assertRaisesRegexp(AssertionError, 'form_insertion_depth.*integer', _convert_input, **polynomials_to_decompose_noninteger_insertion_depth)

    #@attr('active')
    def test_input_check_exponent(self):
        args = self.correct_input.copy()

        args['polynomials_to_decompose'] = ['(a * z0) ** (eps + z1)']
        self.assertRaisesRegexp(AssertionError, 'exponents.*not depend on the .integration_variables', _convert_input, **args)

        args['polynomials_to_decompose'] = ['(a * z0) ** (DummyFunction(eps) + 5)']
        self.assertRaisesRegexp(sp.PolynomialError, 'polynomials.*regulators.*Error while checking: "\( \+ \(a\)\*z0\)\*\*\(DummyFunction\(eps\) \+ 5\)"', _convert_input, **args)

    #@attr('active')
    def test_validate(self):
        self.assertRaisesRegexp(NameError, 'not start with..SecDecInternal', _validate, 'SecDecInternalFunction')
        self.assertRaisesRegexp(NameError, '1a.*cannot be used', _validate, '1a')
        self.assertRaisesRegexp(NameError, 'my_symbol.*cannot be used', _validate, 'my_symbol')
        _validate('symbol1') # should be ok

# --------------------------------- write FORM code ---------------------------------
class TestMakeFORMDefinition(unittest.TestCase):
    #@attr('active')
    def test_function(self):
        polysymbols = sp.symbols("x y z")
        x = Polynomial([[1,0,0]], [1], polysymbols)
        y = Polynomial([[0,1,0]], [1], polysymbols)
        z = Polynomial([[0,0,1]], [1], polysymbols)
        f_dummy = Function('f', x, y, z)
        f = x**2 + y*z

        FORM_code = _make_FORM_definition(f_dummy.symbol, x*x + 3*y*z)
        target_FORM_code = '#define f " + (3)*y*z + (1)*x^2"\n'

        self.assertEqual(FORM_code, target_FORM_code)

#@attr('active')
class TestMakeFORMFunctionDefinition(unittest.TestCase):
    #@attr('active')
    def test_no_args(self):
        symbols = ['x','y']
        x = Polynomial([[1,0]], [1], symbols)
        y = Polynomial([[0,1]], [1], symbols)

        name = 'symbol'
        expression = Sum(x**2, y**2)
        limit = 20
        FORM_code = _make_FORM_function_definition(name, expression, None, limit)

        target_FORM_code  = "  Id symbol = SecDecInternalfDUMMYsymbolPart0+SecDecInternalfDUMMYsymbolPart1;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYsymbolPart0 =  + (1)*x**2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYsymbolPart1 =  + (1)*y**2;\n"

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_sum(self):
        symbols = ['x','y']
        x = Polynomial([[1,0]], [1], symbols)
        y = Polynomial([[0,1]], [1], symbols)

        name = 'myName'
        expression = Sum(Sum(x**2 + 10 * y, x**2 + 10 * y), y * x)
        limit = 20
        FORM_code = _make_FORM_function_definition(name, expression, symbols, limit)

        target_FORM_code  = "  Id myName(x?,y?) = SecDecInternalfDUMMYmyNamePart0(x,y)+SecDecInternalfDUMMYmyNamePart1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x,y)+SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x?,y?) =  + (10)*y + (1)*x**2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x?,y?) =  + (10)*y + (1)*x**2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart1(x?,y?) =  + (1)*x*y;\n"

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_sum_product(self):
        symbols = ['x','y']
        x = Polynomial([[1,0]], [1], symbols)
        y = Polynomial([[0,1]], [1], symbols)

        name = 'myName'
        expression = Sum(Product(x + 2 * y, x**2 + 10 * y), y * x)
        limit = 20
        FORM_code = _make_FORM_function_definition(name, expression, symbols, limit)

        target_FORM_code  = "  Id myName(x?,y?) = SecDecInternalfDUMMYmyNamePart0(x,y)+SecDecInternalfDUMMYmyNamePart1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x,y)*SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x?,y?) =  + (2)*y + (1)*x;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x?,y?) =  + (10)*y + (1)*x**2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart1(x?,y?) =  + (1)*x*y;\n"

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_product_rule(self):
        symbols = ['x','y']
        x = Polynomial([[1,0]], [1], symbols)
        y = Polynomial([[0,1]], [1], symbols)

        name = 'myName'
        expression = ProductRule(Sum(x**2 + 10 * y, y * x), x**2 + 10 * y)
        limit = 20
        FORM_code = _make_FORM_function_definition(name, expression, symbols, limit)

        target_FORM_code  = "  Id myName(x?,y?) = SecDecInternalfDUMMYmyNamePart0(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x,y)*SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part0(x,y)*SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part0Part0(x,y)+SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part0Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part0Part0(x?,y?) =  + (10)*y + (1)*x**2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part0Part1(x?,y?) =  + (1)*x*y;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0Part1(x?,y?) =  + (10)*y + (1)*x**2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x?,y?) =  + (1);\n"

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_polynomial(self):
        symbols = ['x','y']
        x = Polynomial([[1,0]], [1], symbols)
        y = Polynomial([[0,1]], [1], symbols)

        name = 'myName'
        expression = (x**2 + 10 * y + y * x) * (x**2 + 10 * y)
        limit = 20
        FORM_code = _make_FORM_function_definition(name, expression, symbols, limit)

        # ``expression`` has type `Polynomial` --> fall back to rescue since splitting is not implemented
        target_FORM_code  = "  Id myName(x?,y?) =  + (100)*y**2 + (10)*x*y**2 + (20)*x**2*y + (1)*x**3*y + (1)*x**4;\n"

        self.assertEqual(FORM_code, target_FORM_code)

#@attr('active')
class TestMakeFORMSeriesInitialization(unittest.TestCase):
    #@attr('active')
    def test_one_variable(self):
        min_orders = [-2]
        max_orders = [+2]
        sector_ID = 42

        FORM_code = _make_FORM_Series_initilization(min_orders, max_orders, sector_ID, contour_deformation=False)

        target_FORM_code  = '{-2,2,{'
        target_FORM_code +=    '{42,sector_42_order_n2_numIV,sector_42_order_n2_integrand},'
        target_FORM_code +=    '{42,sector_42_order_n1_numIV,sector_42_order_n1_integrand},'
        target_FORM_code +=    '{42,sector_42_order_0_numIV,sector_42_order_0_integrand},'
        target_FORM_code +=    '{42,sector_42_order_1_numIV,sector_42_order_1_integrand},'
        target_FORM_code +=    '{42,sector_42_order_2_numIV,sector_42_order_2_integrand}'
        target_FORM_code += '},true}'

        print('is:')
        print(FORM_code)
        print()
        print('should:')
        print(target_FORM_code)
        print('----------------')
        print()

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_two_variables(self):
        min_orders = [-2, -1]
        max_orders = [+0, +2]
        sector_ID = 8

        FORM_code = _make_FORM_Series_initilization(min_orders, max_orders, sector_ID, contour_deformation=False)

        target_FORM_code  = '{-2,0,{'
        target_FORM_code +=   '{-1,2,{'
        target_FORM_code +=       '{8,sector_8_order_n2_n1_numIV,sector_8_order_n2_n1_integrand},'
        target_FORM_code +=       '{8,sector_8_order_n2_0_numIV,sector_8_order_n2_0_integrand},'
        target_FORM_code +=       '{8,sector_8_order_n2_1_numIV,sector_8_order_n2_1_integrand},'
        target_FORM_code +=       '{8,sector_8_order_n2_2_numIV,sector_8_order_n2_2_integrand}'
        target_FORM_code +=   '},true},'
        target_FORM_code +=   '{-1,2,{'
        target_FORM_code +=       '{8,sector_8_order_n1_n1_numIV,sector_8_order_n1_n1_integrand},'
        target_FORM_code +=       '{8,sector_8_order_n1_0_numIV,sector_8_order_n1_0_integrand},'
        target_FORM_code +=       '{8,sector_8_order_n1_1_numIV,sector_8_order_n1_1_integrand},'
        target_FORM_code +=       '{8,sector_8_order_n1_2_numIV,sector_8_order_n1_2_integrand}'
        target_FORM_code +=   '},true},'
        target_FORM_code +=   '{-1,2,{'
        target_FORM_code +=       '{8,sector_8_order_0_n1_numIV,sector_8_order_0_n1_integrand},'
        target_FORM_code +=       '{8,sector_8_order_0_0_numIV,sector_8_order_0_0_integrand},'
        target_FORM_code +=       '{8,sector_8_order_0_1_numIV,sector_8_order_0_1_integrand},'
        target_FORM_code +=       '{8,sector_8_order_0_2_numIV,sector_8_order_0_2_integrand}'
        target_FORM_code +=   '},true}'
        target_FORM_code += '},true}'

        print('is:')
        print(FORM_code)
        print()
        print('should:')
        print(target_FORM_code)
        print('----------------')
        print()

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_three_variables(self):
        min_orders = [-1, -3, +0]
        max_orders = [+0, -1, +2]
        sector_ID = 90

        FORM_code = _make_FORM_Series_initilization(min_orders, max_orders, sector_ID, contour_deformation=False)

        target_FORM_code  = '{-1,0,{'
        target_FORM_code +=   '{-3,-1,{'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_n1_n3_0_numIV,sector_90_order_n1_n3_0_integrand},'
        target_FORM_code +=       '{90,sector_90_order_n1_n3_1_numIV,sector_90_order_n1_n3_1_integrand},'
        target_FORM_code +=       '{90,sector_90_order_n1_n3_2_numIV,sector_90_order_n1_n3_2_integrand}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_n1_n2_0_numIV,sector_90_order_n1_n2_0_integrand},'
        target_FORM_code +=       '{90,sector_90_order_n1_n2_1_numIV,sector_90_order_n1_n2_1_integrand},'
        target_FORM_code +=       '{90,sector_90_order_n1_n2_2_numIV,sector_90_order_n1_n2_2_integrand}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_n1_n1_0_numIV,sector_90_order_n1_n1_0_integrand},'
        target_FORM_code +=       '{90,sector_90_order_n1_n1_1_numIV,sector_90_order_n1_n1_1_integrand},'
        target_FORM_code +=       '{90,sector_90_order_n1_n1_2_numIV,sector_90_order_n1_n1_2_integrand}'
        target_FORM_code +=     '},true}'
        target_FORM_code +=   '},true},'
        target_FORM_code +=   '{-3,-1,{'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_0_n3_0_numIV,sector_90_order_0_n3_0_integrand},'
        target_FORM_code +=       '{90,sector_90_order_0_n3_1_numIV,sector_90_order_0_n3_1_integrand},'
        target_FORM_code +=       '{90,sector_90_order_0_n3_2_numIV,sector_90_order_0_n3_2_integrand}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_0_n2_0_numIV,sector_90_order_0_n2_0_integrand},'
        target_FORM_code +=       '{90,sector_90_order_0_n2_1_numIV,sector_90_order_0_n2_1_integrand},'
        target_FORM_code +=       '{90,sector_90_order_0_n2_2_numIV,sector_90_order_0_n2_2_integrand}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_0_n1_0_numIV,sector_90_order_0_n1_0_integrand},'
        target_FORM_code +=       '{90,sector_90_order_0_n1_1_numIV,sector_90_order_0_n1_1_integrand},'
        target_FORM_code +=       '{90,sector_90_order_0_n1_2_numIV,sector_90_order_0_n1_2_integrand}'
        target_FORM_code +=     '},true}'
        target_FORM_code +=   '},true}'
        target_FORM_code += '},true}'

        print('is:')
        print(FORM_code)
        print()
        print('should:')
        print(target_FORM_code)
        print('----------------')
        print()

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_three_variables_with_contour_deformation(self):
        min_orders = [-1, -3, +0]
        max_orders = [+0, -1, +2]
        sector_ID = 90

        FORM_code = _make_FORM_Series_initilization(min_orders, max_orders, sector_ID, contour_deformation=True)

        target_FORM_code  = '{-1,0,{'
        target_FORM_code +=   '{-3,-1,{'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_n1_n3_0_numIV,sector_90_order_n1_n3_0_integrand,\n' + \
        '                           sector_90_order_n1_n3_0_contour_deformation,sector_90_order_n1_n3_0_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n3_0_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_n1_n3_1_numIV,sector_90_order_n1_n3_1_integrand,\n' + \
        '                           sector_90_order_n1_n3_1_contour_deformation,sector_90_order_n1_n3_1_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n3_1_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_n1_n3_2_numIV,sector_90_order_n1_n3_2_integrand,\n' + \
        '                           sector_90_order_n1_n3_2_contour_deformation,sector_90_order_n1_n3_2_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n3_2_maximal_allowed_deformation_parameters}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_n1_n2_0_numIV,sector_90_order_n1_n2_0_integrand,\n' + \
        '                           sector_90_order_n1_n2_0_contour_deformation,sector_90_order_n1_n2_0_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n2_0_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_n1_n2_1_numIV,sector_90_order_n1_n2_1_integrand,\n' + \
        '                           sector_90_order_n1_n2_1_contour_deformation,sector_90_order_n1_n2_1_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n2_1_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_n1_n2_2_numIV,sector_90_order_n1_n2_2_integrand,\n' + \
        '                           sector_90_order_n1_n2_2_contour_deformation,sector_90_order_n1_n2_2_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n2_2_maximal_allowed_deformation_parameters}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_n1_n1_0_numIV,sector_90_order_n1_n1_0_integrand,\n' + \
        '                           sector_90_order_n1_n1_0_contour_deformation,sector_90_order_n1_n1_0_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n1_0_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_n1_n1_1_numIV,sector_90_order_n1_n1_1_integrand,\n' + \
        '                           sector_90_order_n1_n1_1_contour_deformation,sector_90_order_n1_n1_1_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n1_1_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_n1_n1_2_numIV,sector_90_order_n1_n1_2_integrand,\n' + \
        '                           sector_90_order_n1_n1_2_contour_deformation,sector_90_order_n1_n1_2_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_n1_n1_2_maximal_allowed_deformation_parameters}'
        target_FORM_code +=     '},true}'
        target_FORM_code +=   '},true},'
        target_FORM_code +=   '{-3,-1,{'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_0_n3_0_numIV,sector_90_order_0_n3_0_integrand,\n' + \
        '                           sector_90_order_0_n3_0_contour_deformation,sector_90_order_0_n3_0_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n3_0_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_0_n3_1_numIV,sector_90_order_0_n3_1_integrand,\n' + \
        '                           sector_90_order_0_n3_1_contour_deformation,sector_90_order_0_n3_1_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n3_1_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_0_n3_2_numIV,sector_90_order_0_n3_2_integrand,\n' + \
        '                           sector_90_order_0_n3_2_contour_deformation,sector_90_order_0_n3_2_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n3_2_maximal_allowed_deformation_parameters}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_0_n2_0_numIV,sector_90_order_0_n2_0_integrand,\n' + \
        '                           sector_90_order_0_n2_0_contour_deformation,sector_90_order_0_n2_0_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n2_0_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_0_n2_1_numIV,sector_90_order_0_n2_1_integrand,\n' + \
        '                           sector_90_order_0_n2_1_contour_deformation,sector_90_order_0_n2_1_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n2_1_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_0_n2_2_numIV,sector_90_order_0_n2_2_integrand,\n' + \
        '                           sector_90_order_0_n2_2_contour_deformation,sector_90_order_0_n2_2_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n2_2_maximal_allowed_deformation_parameters}'
        target_FORM_code +=     '},true},'
        target_FORM_code +=     '{0,2,{'
        target_FORM_code +=       '{90,sector_90_order_0_n1_0_numIV,sector_90_order_0_n1_0_integrand,\n' + \
        '                           sector_90_order_0_n1_0_contour_deformation,sector_90_order_0_n1_0_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n1_0_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_0_n1_1_numIV,sector_90_order_0_n1_1_integrand,\n' + \
        '                           sector_90_order_0_n1_1_contour_deformation,sector_90_order_0_n1_1_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n1_1_maximal_allowed_deformation_parameters},'
        target_FORM_code +=       '{90,sector_90_order_0_n1_2_numIV,sector_90_order_0_n1_2_integrand,\n' + \
        '                           sector_90_order_0_n1_2_contour_deformation,sector_90_order_0_n1_2_contour_deformation_polynomial,\n' + \
        '                           sector_90_order_0_n1_2_maximal_allowed_deformation_parameters}'
        target_FORM_code +=     '},true}'
        target_FORM_code +=   '},true}'
        target_FORM_code += '},true}'

        print('is:')
        print(FORM_code)
        print()
        print('should:')
        print(target_FORM_code)
        print('----------------')
        print()

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_three_variables_simple(self):
        min_orders = [+0, -1, +2]
        max_orders = [+0, -1, +2]
        sector_ID = 90

        FORM_code = _make_FORM_Series_initilization(min_orders, max_orders, sector_ID, contour_deformation=False)

        target_FORM_code  = '{0,0,{'
        target_FORM_code +=   '{-1,-1,{'
        target_FORM_code +=     '{2,2,{'
        target_FORM_code +=         '{90,sector_90_order_0_n1_2_numIV,sector_90_order_0_n1_2_integrand}'
        target_FORM_code +=     '},true}'
        target_FORM_code +=   '},true}'
        target_FORM_code += '},true}'

        print('is:')
        print(FORM_code)
        print()
        print('should:')
        print(target_FORM_code)
        print('----------------')
        print()

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_three_variables_simple_with_contour_deformation(self):
        min_orders = [+0, -1, +2]
        max_orders = [+0, -1, +2]
        sector_ID = 90

        FORM_code = _make_FORM_Series_initilization(min_orders, max_orders, sector_ID, contour_deformation=True)

        target_FORM_code  = '{0,0,{'
        target_FORM_code +=   '{-1,-1,{'
        target_FORM_code +=     '{2,2,{'
        target_FORM_code +=         '{90,sector_90_order_0_n1_2_numIV,sector_90_order_0_n1_2_integrand,\n' + \
          '                           sector_90_order_0_n1_2_contour_deformation,sector_90_order_0_n1_2_contour_deformation_polynomial,\n' + \
          '                           sector_90_order_0_n1_2_maximal_allowed_deformation_parameters}'
        target_FORM_code +=     '},true}'
        target_FORM_code +=   '},true}'
        target_FORM_code += '},true}'

        print('is:')
        print(FORM_code)
        print()
        print('should:')
        print(target_FORM_code)
        print('----------------')
        print()

        self.assertEqual(FORM_code, target_FORM_code)

class TestMiscellaneous(unittest.TestCase):
    #@attr('active')
    def test_derivative_muliindex_to_name(self):
        basename = 'f'
        multiindex = (1,2,1)

        result = _derivative_muliindex_to_name(basename, multiindex)
        target_result = 'ddddfd0d1d1d2'
        self.assertEqual(result, target_result)

    #@attr('active')
    def test_make_FORM_list(self):
        python_list = ['a', 'b', 'c']
        FORM_list = _make_FORM_list(python_list)
        target_FORM_list = 'a,b,c'
        self.assertEqual(FORM_list, target_FORM_list)

    #@attr('active')
    def test_make_FORM_shifted_orders(self):
        powers = [(0,0,0), (1,0,0), (0,1,1)]

        FORM_code = _make_FORM_shifted_orders(powers)

        target_FORM_code  = '#define shiftedRegulator1PowerOrder1 "0"\n'
        target_FORM_code += '#define shiftedRegulator2PowerOrder1 "0"\n'
        target_FORM_code += '#define shiftedRegulator3PowerOrder1 "0"\n'

        target_FORM_code += '#define shiftedRegulator1PowerOrder2 "1"\n'
        target_FORM_code += '#define shiftedRegulator2PowerOrder2 "0"\n'
        target_FORM_code += '#define shiftedRegulator3PowerOrder2 "0"\n'

        target_FORM_code += '#define shiftedRegulator1PowerOrder3 "0"\n'
        target_FORM_code += '#define shiftedRegulator2PowerOrder3 "1"\n'
        target_FORM_code += '#define shiftedRegulator3PowerOrder3 "1"'

        self.assertEqual(FORM_code, target_FORM_code)

class TestWriteCppCodePrefactor(unittest.TestCase):
    #@attr('active')
    def test_one_regulator(self):
        expanded_prefactor = Polynomial([[-1],[0],[1]],['-c0','r0','r1'], ['eps'])
        real_parameters = sp.sympify(['r0','r1'])
        complex_parameters = sp.sympify(['c0'])

        for i in range(2):
            if i == 0:
                expanded_prefactor.truncated = True
            else:
                expanded_prefactor.truncated = False

            cpp_code = _make_prefactor_function(expanded_prefactor, real_parameters, complex_parameters)

            target_cpp_code  =         '#define r0 real_parameters.at(0)\n'
            target_cpp_code += '        #define r1 real_parameters.at(1)\n'
            target_cpp_code += '        #define c0 complex_parameters.at(0)\n'
            if i == 0:
                target_cpp_code += '        return {-1,1,{{-c0},{r0},{r1}},true};\n'
            else:
                target_cpp_code += '        return {-1,1,{{-c0},{r0},{r1}},false};\n'
            target_cpp_code += '        #undef r0\n'
            target_cpp_code += '        #undef r1\n'
            target_cpp_code += '        #undef c0'

            print('i =', i)
            print('cpp_code')
            print(cpp_code)
            print()
            print('target_cpp_code')
            print(target_cpp_code)
            print('-----')

            self.assertEqual(cpp_code, target_cpp_code)

    #@attr('active')
    def test_two_regulators(self):
        symbols = sp.sympify(['alpha','eps'])
        alpha_coeffs = [
                           Polynomial([[0,0],[0,1]], ['r0','c1'], symbols),
                           Polynomial([[0,-1]], ['c1'], symbols),
                           Polynomial([[0,1],[0,2]], ['c0','c1'], symbols)
                       ]
        expanded_prefactor = Polynomial([[-1,0],[0,0],[1,0]], alpha_coeffs, symbols)
        real_parameters = sp.sympify(['r0'])
        complex_parameters = sp.sympify(['c0','c1'])

        true_or_false = lambda b: 'true' if b else 'false'

        for i in range(3):
            if i == 0:
                expanded_prefactor.truncated = False
                for coeff in expanded_prefactor.coeffs:
                    coeff.truncated = True
            elif i == 1:
                expanded_prefactor.truncated = True
                for coeff in expanded_prefactor.coeffs:
                    coeff.truncated = False
            else:
                expanded_prefactor.truncated = True
                expanded_prefactor.coeffs[0].truncated = False
                expanded_prefactor.coeffs[1].truncated = True
                expanded_prefactor.coeffs[2].truncated = False

            cpp_code = _make_prefactor_function(expanded_prefactor, real_parameters, complex_parameters)

            target_cpp_code  =         '#define r0 real_parameters.at(0)\n'
            target_cpp_code += '        #define c0 complex_parameters.at(0)\n'
            target_cpp_code += '        #define c1 complex_parameters.at(1)\n'

            target_cpp_code += '        return {-1,1,{'
            target_cpp_code +=             '{0,1,{{r0},{c1}},%s},' % true_or_false(expanded_prefactor.coeffs[0].truncated)
            target_cpp_code +=             '{-1,-1,{{c1}},%s},' % true_or_false(expanded_prefactor.coeffs[1].truncated)
            target_cpp_code +=             '{1,2,{{c0},{c1}},%s}' % true_or_false(expanded_prefactor.coeffs[2].truncated)
            target_cpp_code +=         '},%s};\n' % true_or_false(expanded_prefactor.truncated)

            target_cpp_code += '        #undef r0\n'
            target_cpp_code += '        #undef c0\n'
            target_cpp_code += '        #undef c1'

            print('i =', i)
            print(cpp_code)
            print()
            print(target_cpp_code)
            print()
            print('-------------------------')

            self.assertEqual(cpp_code, target_cpp_code)

class TestWriteCppCodeFunctionDeclaration(unittest.TestCase):
    #@attr('active')
    def test_zero_args(self):
        code = _make_CXX_function_declaration(function_name = 'f', number_of_arguments = 0)
        target_code = '    integrand_return_t f();\n'
        self.assertEqual(code, target_code)

    #@attr('active')
    def test_one_arg(self):
        code = _make_CXX_function_declaration(function_name = 'f', number_of_arguments = 1)

        target_code  = '    template<typename T0>\n'
        target_code += '    integrand_return_t f(T0 arg0);\n'

        self.assertEqual(code, target_code)

    #@attr('active')
    def test_two_args(self):
        code = _make_CXX_function_declaration(function_name = 'f', number_of_arguments = 2)

        target_code  = '    template<typename T0, typename T1>\n'
        target_code += '    integrand_return_t f(T0 arg0, T1 arg1);\n'

        self.assertEqual(code, target_code)
