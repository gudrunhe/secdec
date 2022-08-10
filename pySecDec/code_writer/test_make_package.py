from .make_package import *
from .make_package import _convert_input, _make_FORM_definition, \
                          _make_FORM_function_definition, _make_FORM_list, \
                          _derivative_muliindex_to_name, _make_FORM_shifted_orders, \
                          _validate, _make_prefactor_function, \
                          _make_CXX_function_declaration
from ..algebra import Function, Polynomial, Product, ProductRule, Sum
from ..misc import sympify_expression
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

# --------------------------------- mid-level tests ---------------------------------
class TestMidLevel(TestMakePackage):
    def check_no_polynomial_names_symmetry(self, split=False, primary_decomposition=False):
        self.tmpdir = 'tmpdir_test_no_polynomial_names_symmetry_python' + python_major_version + '_split_' + str(split) + '_primary_' + str(primary_decomposition)

        # The following input generates a fake-symmetry between ``x``, ``y`` and ``p1``, ``p2``.
        # However, the `polynomial_names` ``p1`` and ``p2`` should be excluded from all symmetry
        # finders.
        template_replacements = \
        make_package(
                        name=self.tmpdir,
                        integration_variables = ['x','y','a','b'],
                        regulators = ['eps','alpha'],

                        requested_orders = [0,0],
                        polynomials_to_decompose = ['(x+y)^(-2+alpha)','x^2+y^2'],
                        polynomial_names=['p1','p2'],
                        other_polynomials=['x*p1+y*p2+a+eps'],

                        decomposition_method='iterative' if primary_decomposition else 'iterative_no_primary',

                        use_Pak=True,
                        use_dreadnaut=True,

                        split=split
                    )

        self.assertEqual(template_replacements['number_of_sectors'],6 if primary_decomposition else 2) # should ignore the symmetry involving `polynomial_names`

    #@attr('active')
    @attr('slow')
    def test_no_polynomial_names_symmetry(self):
        self.check_no_polynomial_names_symmetry()

    #@attr('active')
    @attr('slow')
    def test_no_polynomial_names_symmetry_split(self):
        self.check_no_polynomial_names_symmetry(split=True)

    #@attr('active')
    @attr('slow')
    def test_no_polynomial_names_symmetry_primary(self):
        self.check_no_polynomial_names_symmetry(primary_decomposition=True)

    #@attr('active')
    @attr('slow')
    def test_no_polynomial_names_symmetry_split_primary(self):
        self.check_no_polynomial_names_symmetry(split=True, primary_decomposition=True)

    #@attr('active')
    def test_pole_structures(self):
        self.tmpdir = 'tmpdir_test_pole_structures_python' + python_major_version

        template_replacements = \
        make_package(
                        name=self.tmpdir,
                        integration_variables = ['x','y'],
                        regulators = ['eps'],

                        requested_orders = [0],
                        polynomials_to_decompose = ['(x+y)^(-2+eps)']
                    )

        self.assertEqual(template_replacements['pole_structures_initializer'], '{{-1,0}}')

    #@attr('active')
    def test_pole_structures_serial(self):
        self.tmpdir = 'tmpdir_test_pole_structures_python' + python_major_version

        template_replacements = \
        make_package(
                        name=self.tmpdir,
                        integration_variables = ['x','y'],
                        regulators = ['eps'],
                        processes=1,

                        requested_orders = [0],
                        polynomials_to_decompose = ['(x+y)^(-2+eps)']
                    )

        self.assertEqual(template_replacements['pole_structures_initializer'], '{{-1,0}}')

# ----------------------------------- parse input -----------------------------------
class TestConvertInput(TestMakePackage):
    def setUp(self):
        self.tmpdir = 'tmpdir_test_convert_input_python' + python_major_version
        self.correct_input = dict(
                                      name=self.tmpdir,
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
                                      decomposition_method='iterative_no_primary',
                                      pylink_qmc_transforms=['korobov3x3']
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
        self.assertRaisesRegexp(ValueError, r"\(\-s\) \+ \(\-t\)\*x0\*x1\*x2.*polynomials_to_decompose.*symbols.*\(is.*x0, x1, x2.*should.*z0, z1, z2, eps, alpha", _convert_input, **polynomials_to_decompose_unrelated_polysymbols)

        polynomials_to_decompose_wrong_polysymbols_in_exponent = self.correct_input.copy()
        polynomials_to_decompose_wrong_polysymbols_in_exponent['polynomials_to_decompose'] = ['1', ExponentiatedPolynomial([[0,0,0,0,0,0,0],[1,1,1,0,0,0,0]],['-s','-t'],polysymbols=['z0','z1','z2','eps','alpha','U','F'],exponent=Polynomial([[0,0,1]],[1],['x0','x1','x2']))]
        _convert_input(**polynomials_to_decompose_wrong_polysymbols_in_exponent) # should be ok

        polynomials_to_decompose_polynomial_in_coeff = self.correct_input.copy()
        polynomials_to_decompose_polynomial_in_coeff['polynomials_to_decompose'] = [1,Polynomial([[0,0,0,0,0,0,0],[1,1,1,0,0,0,0]],[Polynomial([[3,2,0,0,0,0,0]], ['-s'], ['z0','z1','z2','eps','alpha','U','F']),'-t'],['z0','z1','z2','eps','alpha','U','F'])]
        _convert_input(**polynomials_to_decompose_polynomial_in_coeff) # should be ok

        polynomials_to_decompose_sympy_exponent = self.correct_input.copy()
        polynomials_to_decompose_sympy_exponent['polynomials_to_decompose'] = ['1', ExponentiatedPolynomial([[0,0,0,0,0,0,0],[1,1,1,0,0,0,0]],['-s','-t'],polysymbols=['z0','z1','z2','eps','alpha','U','F'],exponent='1+eps')]
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
    def test_validate_basic(self):
        self.assertRaisesRegexp(NameError, 'not begin with.*SecDecInternal', _validate, 'SecDecInternalFunction')
        self.assertRaisesRegexp(NameError, '1a.*cannot be used', _validate, '1a')
        self.assertRaisesRegexp(NameError, 'my_symbol.*cannot contain.*underscore.*_', _validate, 'my_symbol')
        _validate('symbol1') # should be ok

    #@attr('active')
    def test_validate_allow_underscore(self):
        self.assertRaisesRegexp(NameError, '^"my_name" cannot contain an underscore character "_"$', _validate, 'my_name')
        _validate('my_name', True) # should be ok
        self.assertRaisesRegexp(NameError, '^"with_underscore" cannot contain an underscore character "_"$', _validate, 'with_underscore', allow_underscore=False)
        _validate('with_underscore', allow_underscore=True) # should be ok

    #@attr('active')
    def test_validate_bans(self):
        for allow_underscore in [True, False]:
            self.assertRaisesRegexp(NameError, '^"double" cannot be used as symbol$', _validate, 'double', allow_underscore)
            self.assertRaisesRegexp(NameError, '^"cubareal" cannot be used as symbol$', _validate, 'cubareal', allow_underscore)
            self.assertRaisesRegexp(NameError, '^"float" cannot be used as symbol$', _validate, 'float', allow_underscore)
            self.assertRaisesRegexp(NameError, '^"sqrt" cannot be used as symbol$', _validate, 'sqrt', allow_underscore)
            self.assertRaisesRegexp(NameError, '^"AtomicThing" cannot be used as symbol \(must not begin with "atomic"\)$', _validate, 'AtomicThing', allow_underscore)

        self.assertRaisesRegexp(NameError, '^"_my_name" cannot be used as symbol \(must not begin with "_"\)$', _validate, '_my_name', allow_underscore=True)
        self.assertRaisesRegexp(NameError, '^"_my_name" cannot contain an underscore character "_"$', _validate, '_my_name', allow_underscore=False)

    #@attr('active')
    def test_remainder_expression_with_polynomial_reference(self):
        # `make_package` should raise an error if the `remainder_expression`
        # refers to any of the `polynomial_names`
        keyword_arguments = self.correct_input.copy()
        keyword_arguments['remainder_expression'] = 'firstPolynomialName'
        keyword_arguments['polynomial_names'] = ['firstPolynomialName']
        keyword_arguments['polynomials_to_decompose'] = ['z0 + z1 + z2']
        self.assertRaisesRegexp(ValueError, r'polynomial_names.*firstPolynomialName.*not.*remainder_expression', _convert_input, **keyword_arguments)

    #@attr('active')
    def test_polynomials_to_decompose_self_reference(self):
        # `make_package` should raise an error if any of the `polynomials_to_decompose`
        # refers to any of the `polynomial_names`
        keyword_arguments = self.correct_input.copy()
        keyword_arguments['polynomial_names'] = ['firstPolynomialName']
        keyword_arguments['polynomials_to_decompose'] = ['z0 + firstPolynomialName']
        self.assertRaisesRegexp(ValueError, r'polynomial_names.*firstPolynomialName.*not.*polynomials_to_decompose', _convert_input, **keyword_arguments)

    #@attr('active')
    def test_positive_polynomials(self):
        keyword_arguments = self.correct_input.copy()
        keyword_arguments['positive_polynomials'] = ['U','missingInPolynomialNames','F']
        self.assertRaisesRegexp(AssertionError, r'missingInPolynomialNames.*positive_polynomials.*not.*polynomial_names', _convert_input, **keyword_arguments)

        keyword_arguments = self.correct_input.copy()
        keyword_arguments['positive_polynomials'] = ['U','not_a + symbol','F']
        self.assertRaisesRegexp(AssertionError, r'All.*positive_polynomials.*symbols', _convert_input, **keyword_arguments)

    #@attr('active')
    def test_regulator_in_polynomial_to_decompose(self):
        keyword_arguments = self.correct_input.copy()
        keyword_arguments['polynomials_to_decompose'] = ['z0 + firstPolynomialName', 'z0*z1 + (eps+1)*z1*z2']
        self.assertRaisesRegexp(AssertionError, r'polynomials_to_decompose.*not depend.*regulators', _convert_input, **keyword_arguments)

    #@attr('active')
    def test_selective_ibp_power_goals_error_message(self):
        keyword_arguments = self.correct_input.copy()
        integration_variables = ['z0','z1','z2']
        keyword_arguments['integration_variables'] = list(integration_variables)
        keyword_arguments['ibp_power_goal'] = [-1,0]
        self.assertRaisesRegexp(AssertionError, r'number of .*ibp_power_goal.* \(2\).* (match|equal) .*number of .*integration_variables.* \(3\)', _convert_input, **keyword_arguments)

    #@attr('active')
    def test_selective_ibp_power_goals(self):
        keyword_arguments = self.correct_input.copy()
        integration_variables = ['z0','z1','z2']
        ibp_power_goal = [-1,0,1]
        keyword_arguments['integration_variables'] = list(integration_variables)
        keyword_arguments['ibp_power_goal'] = list(ibp_power_goal)
        converted_input = _convert_input(**keyword_arguments)
        converted_integration_variables = converted_input[1]
        converted_ibp_power_goal = converted_input[2]
        for i in range(3):
            self.assertEqual(converted_integration_variables[i],sp.symbols(integration_variables[i]))
            self.assertEqual(converted_ibp_power_goal[converted_integration_variables[i]], ibp_power_goal[i])

    #@attr('active')
    def test_old_single_ibp_power_goal(self):
        keyword_arguments = self.correct_input.copy()
        integration_variables = ['z0','z1','z2']
        keyword_arguments['integration_variables'] = list(integration_variables)
        ibp_power_goal = 0
        keyword_arguments['ibp_power_goal'] = ibp_power_goal
        converted_input = _convert_input(**keyword_arguments)
        converted_integration_variables = converted_input[1]
        converted_ibp_power_goal = converted_input[2]
        for i in range(3):
            self.assertEqual(converted_integration_variables[i],sp.symbols(integration_variables[i]))
            self.assertEqual(converted_ibp_power_goal[converted_integration_variables[i]], ibp_power_goal)

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
        FORM_code = ''.join(FORM_code)

        target_FORM_code  = "  Id symbol = SecDecInternalfDUMMYsymbolPart0+SecDecInternalfDUMMYsymbolPart1;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYsymbolPart0 =  + (1)*x^2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYsymbolPart1 =  + (1)*y^2;\n"

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
        FORM_code = ''.join(FORM_code)

        target_FORM_code  = "  Id myName(x?,y?) = SecDecInternalfDUMMYmyNamePart0(x,y)+SecDecInternalfDUMMYmyNamePart1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x,y)+SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x?,y?) =  + (10)*y + (1)*x^2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x?,y?) =  + (10)*y + (1)*x^2;\n"
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
        FORM_code = ''.join(FORM_code)

        target_FORM_code  = "  Id myName(x?,y?) = SecDecInternalfDUMMYmyNamePart0(x,y)+SecDecInternalfDUMMYmyNamePart1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x,y)*SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x?,y?) =  + (2)*y + (1)*x;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x?,y?) =  + (10)*y + (1)*x^2;\n"
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
        FORM_code = ''.join(FORM_code)

        target_FORM_code  = "  Id myName(x?,y?) = SecDecInternalfDUMMYmyNamePart0(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x,y)*SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x,y)*SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part2(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x?,y?) =  + (1);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x?,y?) = SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1Part0(x,y)+SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1Part1(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1Part0(x?,y?) =  + (10)*y + (1)*x^2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1Part1(x?,y?) =  + (1)*x*y;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part2(x?,y?) =  + (10)*y + (1)*x^2;\n"

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
        FORM_code = ''.join(FORM_code)

        target_FORM_code  = "  Id myName(x?,y?) = SecDecInternalfDUMMYmyNamePart0(x,y)+SecDecInternalfDUMMYmyNamePart1(x,y)+SecDecInternalfDUMMYmyNamePart2(x,y)+SecDecInternalfDUMMYmyNamePart3(x,y)+SecDecInternalfDUMMYmyNamePart4(x,y);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?) =  + (100)*y^2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart1(x?,y?) =  + (10)*x*y^2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart2(x?,y?) =  + (20)*x^2*y;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart3(x?,y?) =  + (1)*x^3*y;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart4(x?,y?) =  + (1)*x^4;\n"

        self.assertEqual(FORM_code, target_FORM_code)

    #@attr('active')
    def test_polynomial_coefficient(self):
        symbols = ['x','y','z']

        name = 'myName'
        coeff = Polynomial([[1,1,1],[1,1,2],[1,1,3],[1,1,4]], [1,2,3,4], ['a','b','c'])
        expression = Polynomial([[0,1,1]], [coeff], symbols)
        limit = 20
        FORM_code = _make_FORM_function_definition(name, expression, symbols, limit)
        FORM_code = ''.join(FORM_code)

        target_FORM_code  = "  Id myName(x?,y?,z?) = SecDecInternalfDUMMYmyNamePart0(x,y,z) * ( + (1)*y*z);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYmyNamePart0(x?,y?,z?) = SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x,y,z)+SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x,y,z)+SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part2(x,y,z)+SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part3(x,y,z);\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part0(x?,y?,z?) =  + (1)*a*b*c;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part1(x?,y?,z?) =  + (2)*a*b*c^2;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part2(x?,y?,z?) =  + (3)*a*b*c^3;\n"
        target_FORM_code += "  Id SecDecInternalfDUMMYSecDecInternalfDUMMYmyNamePart0Part3(x?,y?,z?) =  + (4)*a*b*c^4;\n"

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
        real_parameters = sympify_expression(['r0','r1'])
        complex_parameters = sympify_expression(['c0'])
        regulator_names = ['a']

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
                target_cpp_code += '        return {-1,1,{{-c0},{r0},{r1}},true,"eps"};\n'
            else:
                target_cpp_code += '        return {-1,1,{{-c0},{r0},{r1}},false,"eps"};\n'
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
        symbols = sympify_expression(['alpha','eps'])
        alpha_coeffs = [
                           Polynomial([[0,0],[0,1]], ['r0','c1'], symbols),
                           Polynomial([[0,-1]], ['c1'], symbols),
                           Polynomial([[0,1],[0,2]], ['c0','c1'], symbols)
                       ]
        expanded_prefactor = Polynomial([[-1,0],[0,0],[1,0]], alpha_coeffs, symbols)
        real_parameters = sympify_expression(['r0'])
        complex_parameters = sympify_expression(['c0','c1'])

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
            target_cpp_code +=             '{0,1,{{r0},{c1}},%s,"eps"},' % true_or_false(expanded_prefactor.coeffs[0].truncated)
            target_cpp_code +=             '{-1,-1,{{c1}},%s,"eps"},' % true_or_false(expanded_prefactor.coeffs[1].truncated)
            target_cpp_code +=             '{1,2,{{c0},{c1}},%s,"eps"}' % true_or_false(expanded_prefactor.coeffs[2].truncated)
            target_cpp_code +=         '},%s,"alpha"};\n' % true_or_false(expanded_prefactor.truncated)

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

# --------------------------------- algebra helper ----------------------------------
class TestRealPartFunction(unittest.TestCase):
    def setUp(self):
        self.number_of_polysymbols = 3
        self.polysymbols = ['x%i' % i for i in range(self.number_of_polysymbols)]
        self.variables = [Polynomial.from_expression('x%i' % i, self.polysymbols) for i in range(self.number_of_polysymbols)]

    #@attr('active')
    def test_base_function(self):
        Re_x0 = RealPartFunction('Re', self.variables[0])
        self.assertEqual( sympify_expression(Re_x0) , sympify_expression('Re(x0)') )

    #@attr('active')
    def test_derivatives(self):
        Re_x0 = RealPartFunction('Re', self.variables[0]*self.variables[1]*self.variables[1])
        dRe_x0d1 = Re_x0.derive(1)
        self.assertEqual( sympify_expression(dRe_x0d1.derive(0)) , sympify_expression('Re(2*x1)') )

class TestMaxDegreeFunction(unittest.TestCase):
    def setUp(self):
        self.polynomial = Polynomial([[0,1,3],[-1,2,1]],[1,2])
        self.exponentiated_polynomial = ExponentiatedPolynomial(self.polynomial.expolist, self.polynomial.coeffs, exponent='exponent')
        self.target_maxdegrees = [np.inf,2,3]


        polysymbols = self.polynomial.polysymbols
        self.variables = [Polynomial.from_expression(x, polysymbols) for x in polysymbols]
        self.variables[-1] = self.variables[-1]**2

    #@attr('active')
    def test_derivative_and_copy(self):
        maxdegrees = np.array([0,1,2])
        f = MaxDegreeFunction('f', *self.variables, maxdegrees=maxdegrees)

        self.assertEqual( sympify_expression(f) , sympify_expression('f(x0,x1,x2**2)') )

        derivative = f.derive(0)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('0') )

        derivative = f.derive(1)
        self.assertEqual( sympify_expression(f.derive(1)) , sympify_expression('dfd1(x0,x1,x2**2)') )
        derivative = derivative.derive(1)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('0') )

        derivative = f.derive(2)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('2*x2*dfd2(x0,x1,x2**2)') )
        derivative = derivative.derive(2)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('2*dfd2(x0,x1,x2**2)+4*x2**2*ddfd2d2(x0,x1,x2**2)') )
        derivative = derivative.derive(2)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('4*x2*ddfd2d2(x0,x1,x2**2)+8*x2*ddfd2d2(x0,x1,x2**2)') )
        derivative = derivative.derive(2).copy()
        self.assertEqual( sympify_expression(derivative) , sympify_expression('4*ddfd2d2(x0,x1,x2**2)+8*ddfd2d2(x0,x1,x2**2)') )
        derivative = derivative.derive(2)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('0') )

        derivative = f.derive(2).derive(1).copy()
        self.assertEqual( sympify_expression(derivative) , sympify_expression('2*x2*ddfd1d2(x0,x1,x2**2)') )
        derivative = derivative.derive(2)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('2*ddfd1d2(x0,x1,x2**2)+4*x2**2*dddfd1d2d2(x0,x1,x2**2)') )
        self.assertEqual( sympify_expression(derivative.derive(0)) , sympify_expression('0') )
        self.assertEqual( sympify_expression(derivative.derive(1)) , sympify_expression('0') )
        derivative = derivative.derive(2).copy()
        self.assertEqual( sympify_expression(derivative) , sympify_expression('4*x2*dddfd1d2d2(x0,x1,x2**2)+8*x2*dddfd1d2d2(x0,x1,x2**2)') )
        self.assertEqual( sympify_expression(derivative.derive(0)) , sympify_expression('0') )
        self.assertEqual( sympify_expression(derivative.derive(1)) , sympify_expression('0') )
        derivative = derivative.derive(2)
        self.assertEqual( sympify_expression(derivative) , sympify_expression('12*dddfd1d2d2(x0,x1,x2**2)') )
        self.assertEqual( sympify_expression(derivative.derive(0)) , sympify_expression('0') )
        self.assertEqual( sympify_expression(derivative.derive(1)) , sympify_expression('0') )
        self.assertEqual( sympify_expression(derivative.derive(2)) , sympify_expression('0') )

    #@attr('active')
    def test_replace(self):
        maxdegrees = np.array([0,1,2])
        f0 = MaxDegreeFunction('f', *self.variables, maxdegrees=maxdegrees).replace(1,1,remove=True)
        f1 = MaxDegreeFunction('f', *self.variables, maxdegrees=maxdegrees).replace(1,1,remove=False)

        for f in [f0,f1]:
            self.assertEqual( sympify_expression(f) , sympify_expression('f(x0,1,x2**2)') )

        for derivative in [f0.derive(0),f1.derive(0)]:
            self.assertEqual( sympify_expression(derivative) , sympify_expression('0') )

        for f in [f0,f1]:
            derivative = f.derive(-1)
            self.assertEqual( sympify_expression(derivative) , sympify_expression('2*x2*dfd2(x0,1,x2**2)') )
            derivative = derivative.derive(-1).copy()
            self.assertEqual( sympify_expression(derivative) , sympify_expression('2*dfd2(x0,1,x2**2)+4*x2**2*ddfd2d2(x0,1,x2**2)') )
            derivative = derivative.derive(-1).copy()
            self.assertEqual( sympify_expression(derivative) , sympify_expression('4*x2*ddfd2d2(x0,1,x2**2)+8*x2*ddfd2d2(x0,1,x2**2)') )
            derivative = derivative.derive(-1).copy()
            self.assertEqual( sympify_expression(derivative) , sympify_expression('4*ddfd2d2(x0,1,x2**2)+8*ddfd2d2(x0,1,x2**2)') )
            derivative = derivative.derive(-1).copy()
            self.assertEqual( sympify_expression(derivative) , sympify_expression('0') )

    #@attr('active')
    def test_get_maxdegrees(self):
        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.polynomial, ignore_subclass=False),
            self.target_maxdegrees
        )
        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.polynomial, ignore_subclass=True),
            self.target_maxdegrees
        )

        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.polynomial, ignore_subclass=False, indices=[0,2]),
            (self.target_maxdegrees[0],np.inf,self.target_maxdegrees[2])
        )
        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.polynomial, ignore_subclass=True, indices=[0,2]),
            (self.target_maxdegrees[0],np.inf,self.target_maxdegrees[2])
        )

        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.polynomial, ignore_subclass=False, indices=[1]),
            (np.inf,self.target_maxdegrees[1],np.inf)
        )
        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.polynomial, ignore_subclass=True, indices=[1]),
            (np.inf,self.target_maxdegrees[1],np.inf)
        )


        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.exponentiated_polynomial, ignore_subclass=False),
            [np.inf] * 3
        )
        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.exponentiated_polynomial, ignore_subclass=True),
            self.target_maxdegrees
        )

        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.exponentiated_polynomial, ignore_subclass=False, indices=[0,2]),
            [np.inf] * 3
        )
        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.exponentiated_polynomial, ignore_subclass=True, indices=[0,2]),
            (self.target_maxdegrees[0],np.inf,self.target_maxdegrees[2])
        )

        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.exponentiated_polynomial, ignore_subclass=False, indices=[1]),
            [np.inf] * 3
        )
        np.testing.assert_array_equal(
            MaxDegreeFunction.get_maxdegrees(self.exponentiated_polynomial, ignore_subclass=True, indices=[1]),
            (np.inf,self.target_maxdegrees[1],np.inf)
        )
