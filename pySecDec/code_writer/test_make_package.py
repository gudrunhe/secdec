"""Unit tests for the template parser module"""

from .make_package import *
from .make_package import _convert_input
from ..algebra import Polynomial
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
