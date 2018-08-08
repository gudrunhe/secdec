from ..metadata import version, git_id
from .template_parser import parse_template_file, parse_template_tree
from ..algebra import Polynomial, ExponentiatedPolynomial
from ..misc import cached_property
from ..loop_integral.loop_package import loop_package

from time import strftime
from itertools import repeat
import os, sys, shutil
import numpy as np
import sympy as sp

class Coefficient(object):
    def __init__(self, numerator, denominator, order):
        assert (isinstance(numerator, list) == True), "numerator has to be a list of the individual terms"
        assert (isinstance(denominator, list) == True), "denominator has to be a list of the individual terms"
        assert (isinstance(order, int) == True), "order has to be an integer"
        self.numerator = numerator
        self.denominator = denominator
        self.order = order

def sum_package(integral_name, package_generators, generators_args, requested_orders, real_parameters=[],
                complex_parameters=[], coefficients=None):
    r'''
    Decompose, subtract and expand an expression.
    Return it as c++ package.

    :param integral_name:
        string;
        The name of the c++ namepace and the output
        directory.

    :param package_generators:
        iterable;
        List that has package_generators
        (:func:`pySecDec.loop_integral.loop_package`
        or :func:`pySecDec.code_writer.make_package`)
        specified in the order in which they should
        be called on the list of package args.

    :param generators_args:
        iterable;
        dictionaries specifying the arguments for the
        package_generatorrs. See
        :func:`pySecDec.code_writer.make_package` and
        :func:`pySecDec.loop_integral.loop_package` for the
        arguments.

    :param requested_orders:
        iterable;
        Requested orders in the regulators.

    :param real_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as real variables.

    :param complex_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as complex variables.

    :param coefficients:
        list of :class:`pySecDec.code_writer.sum_package.Coefficients`;
        coefficients of the integrals.

    '''
    print('running "sum_package" for ' + integral_name)

    assert len(package_generators) == len(generators_args), \
               "package_generators and generators_args have to have the same length,\
                len(package_generators) = " + len(package_generators) + \
                "len(generators_args) = " + len(generators_args)

    assert (coefficients == None or len(coefficients) == len(generators_args)), \
            "coefficients has to either be None in which case 1 is taken as the \
            default coefficients or coefficients has to have the same length as \
            generator_args"

    # construct the c++ type "nested_series_t"
    # for two regulators, the resulting code should read:
    # "secdecutil::Series<secdecutil::Series<T>>"
    nested_series_type = 'secdecutil::Series<' * len(requested_orders) + 'T' + '>' * len(requested_orders)

    names_makefile = []
    includes = ''
    computes = ''
    for package_generator,generator_args in zip(package_generators,generators_args):
        # define the replacements in the sum_package folder
        name = generator_args["name"]
        names_makefile.append(name)
        includes += '#include "' + name + '/' + name + '.hpp"\n'
        computes += '    COMPUTE(' + name + '_res, ' + name + ', 1., divonne);\n'
    names_makefile = ' '.join(names_makefile)

    sum_of_integrals = '    auto result = \n    '
    if not coefficients:
        sum_of_integrals +=  ' + '.join([ generator_args["name"] + '_res\n    ' for generator_args in generators_args])
    else:
        for coefficient,generator_args in zip(coefficients,generators_args):
            name = generator_args["name"]
            for numerator in coefficient.numerators:
                try:
                    sum_of_integrals += ' + pow(' + str(Polynomial(numerator.expolist,numerator.coeffs,numerator.symbols)) + ',' + sp.printing.ccode(numerator.exponent.evalf(20)) + ')'
                except AttributeError:
                    sum_of_integrals += ' + ' + str(numerator)
            for denominator in coefficient.denominators:
                try:
                    sum_of_integrals += '/ pow(' + str(Polynomial(denominator.expolist,denominator.coeffs,denominator.symbols)) + ',' + sp.printing.ccode(denominator.exponent.evalf(20)) + ')'
                except AttributeError:
                    sum_of_integrals += ' / ' + str(denominator)
            sum_of_integrals += '* ' + name + '_res\n    '
    sum_of_integrals += '    ;'

    parameter_define = ''
    for i,real_parameter in enumerate(real_parameters):
        parameter_define += '    #define ' + real_parameter + ' real_parameters[' + str(i) + ']\n'

    for i,complex_parameter in enumerate(complex_parameters):
        parameter_define += '    #define ' + complex_parameter + ' complex_parameters[' + str(i) + ']\n'

    parameter_undef = ''
    for i,real_parameter in enumerate(real_parameters):
        parameter_undef += '    #undef ' + real_parameter + '\n'

    for i,complex_parameter in enumerate(complex_parameters):
        parameter_undef += '    #undef ' + complex_parameter + '\n'

    replacements_in_files = {
                            'name' : integral_name,
                            'includes' : includes,
                            'computes' : computes,
                            'parameter_define' : parameter_define,
                            'parameter_undef' : parameter_undef,
                            'sum_of_integrals' : sum_of_integrals,
                            'integral_names' : names_makefile,
                            'pySecDec_version' : version,
                            'python_version' : sys.version,
                            'pySecDec_git_id' : git_id,
                            'date_time' : strftime("%a %d %b %Y %H:%M"),
                            'nested_series_type' : nested_series_type
                            }
    filesystem_replacements = { 'integrate_name.cpp' : 'integrate_' + integral_name + '.cpp' }

    # get path to the directory with the template files (path relative to directory with this file: "./templates/")
    template_sources = os.path.join(os.path.split(os.path.abspath(__file__))[0],'templates','sum_package')
    parse_template_tree(template_sources, integral_name, replacements_in_files, filesystem_replacements)

    original_working_directory = os.getcwd()

    try:
        os.chdir(integral_name)
        # call package integrator for every integral

        if coefficients is None:
            coefficients = repeat(None)

        for package_generator, generator_args, coefficient in zip(package_generators, generators_args, coefficients):
            generator_args['real_parameters'] = real_parameters
            generator_args['complex_parameters'] = complex_parameters
            if package_generator == loop_package:
                if coefficient is None:
                    generator_args['requested_order'] = requested_orders[0]
                else:
                    generator_args['requested_order'] = requested_orders[0] - coefficient.orders[0]
            else:
                if coefficient is None:
                    generator_args['requested_orders'] = np.asarray(requested_orders)
                else:
                    generator_args['requested_orders'] = np.asarray(requested_orders) - np.asarray(coefficient.orders)
            package_generator(**generator_args)

    finally:
        os.chdir(original_working_directory)
