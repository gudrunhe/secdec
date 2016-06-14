"""Unit tests for the template parser module"""

from .make_package import *
from .make_package import _convert_input, _make_FORM_definition, \
                          _make_FORM_list, _derivative_muliindex_to_name, \
                          _make_FORM_shifted_orders
from ..algebra import Polynomial, Function
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

    #@attr('active')
    def test_convert_input(self):
        correct_input = dict(
                                target_directory=self.tmpdir,
                                name='some_integral',
                                integration_variables=['z0','z1','z2'],
                                regulators=['eps','alpha'],
                                requested_orders=[1,2],
                                polynomials_to_decompose=[1,Polynomial([[0,0,0],[1,1,1]],['-s','-t'],['z0','z1','z2','eps','alpha'])],
                                polynomial_names=['U','F'],
                                other_polynomials=['U*z1 + F'],
                                prefactor=1,
                                remainder_expression='DummyFunction(z0,eps)',
                                functions=['DummyFunction'],
                                other_variables=['s','t'],
                                form_optimization_level=2,
                                form_work_space='500M',
                                stabilize=False,
                                contour_deformation_polynomial=None,
                                decomposition_method='iterative_no_primary'
                            )

        _convert_input(**correct_input) # should be ok

        requested_orders_wrong_shape = correct_input.copy()
        requested_orders_wrong_shape['requested_orders'] = [[1,1],[0,0]]
        self.assertRaisesRegexp(AssertionError, r'requested_orders.*wrong shape.*is \(2, 2\).*should be \(2,\)', _convert_input, **requested_orders_wrong_shape)

        requested_orders_wrong_length = correct_input.copy()
        requested_orders_wrong_length['requested_orders'] = [1,1,1]
        self.assertRaisesRegexp(AssertionError, 'length.*requested_orders.*match.*length.*regulators', _convert_input, **requested_orders_wrong_length)

        polynomials_to_decompose_unrelated_polysymbols = correct_input.copy()
        polynomials_to_decompose_unrelated_polysymbols['polynomials_to_decompose'] = ['1', Polynomial([[0,0,0],[1,1,1]],['-s','-t'],['x0','x1','x2'])]
        self.assertRaisesRegexp(ValueError, r"\(\-s\) \+ \(\-t\)\*x0\*x1\*x2.*polynomials_to_decompose.*symbols.*\(is.*x0, x1, x2.*should.*\z0, z1, z2, eps, alpha", _convert_input, **polynomials_to_decompose_unrelated_polysymbols)

        polynomials_to_decompose_wrong_polysymbols_in_exponent = correct_input.copy()
        polynomials_to_decompose_wrong_polysymbols_in_exponent['polynomials_to_decompose'] = ['1', ExponentiatedPolynomial([[0,0,0,0,0],[1,1,1,0,0]],['-s','-t'],polysymbols=['z0','z1','z2','eps','alpha'],exponent=Polynomial([[0,0,1]],[1],['x0','x1','x2']))]
        self.assertRaisesRegexp(ValueError, r"exponent.*\(\-s\) \+ \(\-t\)\*z0\*z1\*z2.*polynomials_to_decompose.*symbols.*\(is.*z0, z1, z2.*should.*\z0, z1, z2, eps, alpha", _convert_input, **polynomials_to_decompose_wrong_polysymbols_in_exponent)

        polynomials_to_decompose_sympy_exponent = correct_input.copy()
        polynomials_to_decompose_sympy_exponent['polynomials_to_decompose'] = ['1', ExponentiatedPolynomial([[0,0,0,0,0],[1,1,1,0,0]],['-s','-t'],polysymbols=['z0','z1','z2','eps','alpha'],exponent='1+eps')]
        _convert_input(**polynomials_to_decompose_sympy_exponent) # should be ok

        polynomials_to_decompose_wrong_type_in_exponent = correct_input.copy()
        polynomials_to_decompose_wrong_type_in_exponent['polynomials_to_decompose'] = ['1', ExponentiatedPolynomial([[0,0,0,0,0],[1,1,1,0,0]],['-s','-t'],polysymbols=['z0','z1','z2','eps','alpha'],exponent=ExponentiatedPolynomial([[0,0,1]],[1],polysymbols=['z0','z1','z2','eps','alpha']))]
        self.assertRaisesRegexp(ValueError, r"exponent.*\(\-s\) \+ \(\-t\)\*z0\*z1\*z2.*polynomials_to_decompose.*type.*Polynomial.*ExponentiatedPolynomial", _convert_input, **polynomials_to_decompose_wrong_type_in_exponent)

        polynomials_to_decompose_nontrivial_as_string = correct_input.copy()
        polynomials_to_decompose_nontrivial_as_string['polynomials_to_decompose'] = ['1', '(-s -t*z0*z1*z2)**(2-4*eps+alpha)']
        _convert_input(**polynomials_to_decompose_nontrivial_as_string) # should be ok

# --------------------------------- write FORM code ---------------------------------
#@attr('active')
class TestWriteFORMCode(TestMakePackage):
    def setUp(self):
        self.tmpdir = 'tmpdir_test_convert_input_python' + python_major_version

    #@attr('active')
    def test_make_FORM_definition(self):
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
        target_FORM_list = 'a, b, c'
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
