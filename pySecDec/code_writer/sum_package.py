from ..metadata import version, git_id
from .template_parser import validate_pylink_qmc_transforms, generate_pylink_qmc_macro_dict, parse_template_file, parse_template_tree
from ..misc import sympify_symbols, make_cpp_list, chunks
from .make_package import make_package
import pySecDecContrib

from multiprocessing import Pool
from time import strftime
import json
import numpy as np
import os
import re
import shutil
import subprocess
import sympy as sp
import sys
import tempfile

class Coefficient(object):
    r'''
    Store a coefficient expressed as a product of terms
    in the numerator and a product of terms in the
    denominator.

    :param numerators:
        str, or iterable of str;
        The numerator, or the terms in the numerator.

    :param denominators:
        str, or iterable of str;
        The denominator, or the terms in the denominator.

    :param parameters:
        iterable of strings or sympy symbols;
        The symbols other parameters.

        '''
    def __init__(self, numerators, denominators=(), parameters=()):
        if not np.iterable(parameters):
            raise ValueError("parameters must be iterable")
        if not isinstance(numerators, str):
            numerators = "*".join(f"({f})" for f in numerators) or "1"
        if not isinstance(denominators, str):
            denominators = "*".join(f"({f})" for f in denominators) or ""
        self.expression = f"({numerators})/({denominators})" if denominators else numerators
        self.parameters = parameters

    def leading_orders(self, regulators, maxorder=999999):
        r'''
        Calculate and return the lowest orders of the coefficients
        in the regulators. If the whole coefficient is zero, return
        a list of the maximal orders.

        :param regulators:
            iterable of strings or sympy symbols;
            The symbols denoting the regulators.

        :param maxorder:
            int;
            The maximum order of the regulators; if the actual
            order is higher, this value is used instead.
        '''
        regulators = sympify_symbols(regulators, 'All `regulators` must be symbols.')
        ginsh = os.path.join(pySecDecContrib.dirname, 'bin', 'ginsh')
        # Because ginsh only expands *up to a fixed order*, we
        # need to iteratively increase the order until the
        # expansion is non-zero.
        orders = [0 for reg in regulators]
        rvars = [sp.randprime(2**30, 2**31) for x in self.parameters]
        rregs = [sp.randprime(2**30, 2**31) for reg in regulators]
        while True:
            with tempfile.NamedTemporaryFile("w", encoding="utf8") as f:
                f.write("__START;\n")
                for x, rand in zip(self.parameters, rvars):
                    f.write(f"{x}={rand}:\n")
                f.write("__EXPR=(")
                f.write(self.expression)
                f.write("):\n")
                # Is the whole expression zero?
                f.write("is(subs(__EXPR,{");
                f.write(",".join(f"{reg}=={rand}" for reg, rand in zip(regulators, rregs)))
                f.write("})==0);\n")
                for i, reg in enumerate(regulators):
                    f.write("__EXPRi=subs(__EXPR,{");
                    f.write(",".join(
                            f"{reg}=={rregs[j]}"
                            for j, reg in enumerate(regulators)
                            if j != i
                        ))
                    f.write("}):\n")
                    f.write(f"__EXPRi=series_to_poly(series(__EXPRi,{reg},{orders[i]+1})):\n")
                    # Is the expression zero after expanding in the regulator #i?
                    f.write(f"is(__EXPRi==0);\n")
                    # The leading order in the regulator #i.
                    f.write(f"ldegree(__EXPRi,{reg});\n")
                f.write("exit\n")
                f.flush()
                result = subprocess.check_output([ginsh, f.name], encoding="utf8")
            result = re.sub(r".*__START\n", "", result, flags=re.DOTALL)
            assert "__" not in result
            result = result.splitlines()
            assert len(result) == 1 + len(regulators)*2
            result = [int(line) for line in result]
            if result[0]:
                return [maxorder for reg in regulators]
            done = True
            for i in range(len(regulators)):
                if result[1 + i*2 + 0] and orders[i] <= maxorder:
                    orders[i] += 1
                    done = False
            if done:
                return [result[1 + i*2 + 1] for i in range(len(regulators))]

    def write(self, file):
        """
        Write the coefficient into a given file object.
        """
        file.write(self.expression)

    @staticmethod
    def from_string(expression, exclude_parameters=()):
        expression = expression.replace("**", "^")
        parameters = sorted(list(set(re.findall("[a-zA-Z_][a-zA-Z_0-9]*", expression)) - set(exclude_parameters)))
        return Coefficient(expression, parameters=parameters)

def _generate_one_term(gen_index, sums, complex_parameters, name, package_generator, pylink_qmc_transforms, real_parameters, regulators, replacements_in_files, requested_orders, template_sources):

    sub_name = package_generator.name
    replacements_in_files['sub_integral_name'] = sub_name

    package_generator = package_generator._replace(real_parameters = real_parameters,
                                                   complex_parameters = complex_parameters,
                                                   pylink_qmc_transforms = pylink_qmc_transforms)

    # process coefficients
    lowest_coefficient_orders = {}
    coeffsdir = name + '_data'
    os.makedirs(coeffsdir, exist_ok=True)
    for sumidx, (sum_name, amp) in enumerate(sums.items()):
        if gen_index not in amp: continue
        coefficient = amp[gen_index]
        lowest_orders = coefficient.leading_orders(regulators, maxorder=999999)
        if any(lo < 999999 for lo in lowest_orders):
            lowest_coefficient_orders[sumidx] = lowest_orders
            filename = os.path.join(coeffsdir, f"{sub_name}_coefficient{sumidx}.txt")
            filename2 = os.path.join("disteval", "coefficients", f"{sub_name}_coefficient{sumidx}.txt")
            with open(filename,'w') as coeffile:
                coefficient.write(coeffile)
            os.link(filename, filename2)
    replacements_in_files['lowest_coefficient_orders'] = '{' + '},{'.join(
            ','.join(map(str,lowest_coefficient_orders.get(i, "999999")))
            for i in range(len(sums))
        ) + '}'

    minimal_lowest_coefficient_orders = min(lowest_coefficient_orders.values(), default=999999)
    package_generator = package_generator._replace(regulators = regulators)
    package_generator = package_generator._replace(requested_orders = np.asarray(requested_orders) - minimal_lowest_coefficient_orders)

    # parse integral specific files
    for ch in 'ch':
        suffix = '_weighted_integral.%spp' % ch
        parse_template_file(os.path.join(template_sources, 'src', 'name' + suffix), # source
                            os.path.join('src', sub_name + suffix), # dest
                            replacements_in_files)

    # Call make_package with its arguments
    mp = make_package(**package_generator._asdict())
    highest_coefficient_orders = np.asarray(requested_orders) - mp["description"]["prefactor_lowest_orders"] - mp["description"]["lowest_orders"]
    return lowest_coefficient_orders, highest_coefficient_orders, mp


def sum_package(name, package_generators, regulators, requested_orders,
                real_parameters=[], complex_parameters=[], coefficients=None,
                form_executable=None, pylink_qmc_transforms=['korobov3x3'],
                processes=1):
    r'''
    Decompose, subtract and expand a list of expressions
    of the form

    .. math::
        \sum_j c_{ij} \ \int f_j

    Generate a c++ package with an optmized algorithm to
    evaluate the integrals numerically. It writes the names
    of the integrals in the file `"integral_names.txt"`.
    For the format of the file see :class:`~pySecDec.amplitude_interface.Parser`.

    :param name:
        string;
        The name of the c++ namepace and the output
        directory.

    :param package_generators:
        iterable of `pySecDec.code_writer.MakePackage` and/or
        `pySecDec.loop_integral.LoopPackage` `namedtuples`;
        The generator functions for the integrals
        :math:`\int f_i`
        
        .. note::
            The `pySecDec.code_writer.MakePackage` and
            `pySecDec.loop_integral.LoopPackage` objects
            have the same argument list as their respective
            parent functions
            :func:`pySecDec.code_writer.make_package`
            and :func:`pySecDec.loop_integral.loop_package`.

    :param regulators:
        iterable of strings or sympy symbols;
        The UV/IR regulators of the integral.

    :param requested_orders:
        iterable of integers;
        Compute the expansion in the `regulators` to
        these orders.

    :param real_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as real variables.

    :param complex_parameters:
        iterable of strings or sympy symbols, optional;
        Symbols to be interpreted as complex variables.

    :param coefficients:
        dict of str -> list of strings, optional;
        A map from the (unique) names of the sums to coefficient
        vectors :math:`c_{i}` defining them.

        For backward compatibility, iterable of iterable of
        :class:`.Coefficient`, and iterable of dict of int ->
        :class:`.Coefficient` are accepted as well, with names
        of the sums names set to default values.

    :param pylink_qmc_transforms:
        list or None, optional;
        Required qmc integral transforms, options are:

        * ``none``
        * ``baker``
        * ``korobov<i>x<j>`` for 1 <= i,j <= 6
        * ``korobov<i>`` for 1 <= i <= 6 (same as ``korobov<i>x<i>``)
        * ``sidi<i>`` for 1 <= i <= 6

        Default: ``['korobov3x3']``
    
    :param processes:
        integer or None, optional;
        Parallelize the generation of terms in a sum using this
        many processes.

        When set to a value larger than 1, this will override the
        ``processes`` argument of the terms in a sum, meaning
        that each term will not use parallelization, but rather
        different terms will be generated in parallel.

        Default: ``1``

    '''
    print('running "sum_package" for ' + name)

    # convert input iterables to lists
    package_generators = list(package_generators)
    regulators = list(regulators)
    requested_orders = list(requested_orders)
    real_parameters = list(real_parameters)
    complex_parameters = list(complex_parameters)

    # check that the names of all integrals and the sum itself
    # do not clash
    assert len(set([package_generator.name for package_generator in package_generators] + [name])) == len(package_generators) + 1, \
        "The `name` parameter, and the names of all integrals must all be distinct."

    pylink_qmc_transforms = validate_pylink_qmc_transforms(pylink_qmc_transforms)

    # prepare coefficients
    coefficient_1 = Coefficient((), (), ())
    if coefficients is None:
        sums = {
            name: {i : coefficient_1 for i in range(len(package_generators))}
        }
    elif isinstance(coefficients, dict):
        for sum_name, c_i in coefficients.items():
            if isinstance(c_i, dict):
                for genidx in c_i.keys():
                    if not (0 <= genidx < len(package_generators)):
                        raise ValueError(f"coefficient for generator {genidx} is provided, but there are only {len(package_generators)}")
            else:
                if len(c_i) != len(package_generators):
                    raise ValueError(f"`coefficients` must have list of the same length as `package_generators`.")
        sums = {
            sum_name : {
                genidx : (
                    coefficient_1 if c is None else
                    c if isinstance(c, Coefficient) else
                    Coefficient.from_string(str(c), exclude_parameters=regulators)
                )
                for genidx, c in (c_i.items() if isinstance(c_i, dict) else enumerate(c_i))
                if c != 0 and c != "0"
            }
            for sum_name, c_i in coefficients.items()
        }
    else:
        for i, c_i in enumerate(coefficients):
            if isinstance(c_i, dict):
                for genidx in c_i.keys():
                    if not (0 <= genidx < len(package_generators)):
                        raise ValueError(f"coefficient for generator {genidx} is provided, but there are only {len(package_generators)}")
            else:
                if len(c_i) != len(package_generators):
                    raise ValueError(f"`coefficients` must have list of the same length as `package_generators`.")
        sums = {
            f"sum{sumidx}": {
                genidx : (
                    coefficient_1 if c is None else
                    c if isinstance(c, Coefficient) else
                    Coefficient.from_string(str(c), exclude_parameters=regulators)
                )
                for genidx, c in (c_i.items() if isinstance(c_i, dict) else enumerate(c_i))
                if c != 0 and c != "0"
            }
            for sumidx, c_i in enumerate(coefficients)
        }

    # construct the c++ type "nested_series_t"
    # for two regulators, the resulting code should read:
    # "secdecutil::Series<secdecutil::Series<T>>"
    nested_series_type = 'secdecutil::Series<' * len(requested_orders) + 'T' + '>' * len(requested_orders)

    # no need to check this, package_generator.regulators is ignored
    #for package_generator in package_generators:
    #    assert len(requested_orders) == len(package_generator.regulators), \
    #        "The `requested_orders` must match the number of regulators"

    # define required listings of contributing integrals
    sub_integral_names = []
    integral_initialization = []
    integral_initialization_with_contour_deformation = []
    weighted_integral_includes = []
    weighted_integral_sum_initialization = []
    for package_generator in package_generators:
        sub_name = package_generator.name
        sub_integral_names.append(sub_name)
        integral_initialization.append( 'std::vector<nested_series_t<sum_t>> integral_' + sub_name + ' = ' + sub_name + '::make_integral(real_parameters,complex_parameters,integrator);' )
        integral_initialization_with_contour_deformation.append( 'std::vector<nested_series_t<sum_t>> integral_' + sub_name + ' = ' + sub_name + '::make_integral(real_parameters,complex_parameters,integrator,number_of_presamples,deformation_parameters_maximum,deformation_parameters_minimum,deformation_parameters_decrease_factor);' )
        weighted_integral_includes.append( '#include "' + sub_name + '_weighted_integral.hpp"')
        weighted_integral_sum_initialization.append( 'amplitude += ' + sub_name + '::make_weighted_integral(real_parameters, complex_parameters, integral_' + sub_name + ', amp_idx, lib_path);' )
    sub_integral_names = ' '.join(sub_integral_names)
    integral_initialization = '\n        '.join(integral_initialization)
    integral_initialization_with_contour_deformation = '\n        '.join(integral_initialization_with_contour_deformation)
    weighted_integral_includes = '\n'.join(weighted_integral_includes)
    weighted_integral_sum_initialization[0] = weighted_integral_sum_initialization[0].replace('+',' ',1)
    weighted_integral_sum_initialization = '\n            '.join(weighted_integral_sum_initialization)

    # generate translation from transform short names 'korobov#x#' and 'sidi#' to C++ macros
    pylink_qmc_dynamic_cast_integrator = generate_pylink_qmc_macro_dict('DYNAMIC_CAST_INTEGRATOR')
    pylink_qmc_instantiate_make_amplitudes = generate_pylink_qmc_macro_dict('INSTANTIATE_MAKE_AMPLITUDES')
    pylink_qmc_instantiate_make_integral = generate_pylink_qmc_macro_dict('INSTANTIATE_MAKE_INTEGRAL')
    pylink_qmc_instantiate_amplitude_integral = generate_pylink_qmc_macro_dict('INSTANTIATE_AMPLITUDE_INTEGRAL')
    pylink_qmc_instantiations_translation = generate_pylink_qmc_macro_dict('INSTANTIATE')
    pylink_qmc_extern_translation = generate_pylink_qmc_macro_dict('EXTERN')
    pylink_qmc_case_translation = generate_pylink_qmc_macro_dict('CASE')

    # parse the required pylink templates and generate a list of files to write
    pylink_qmc_dynamic_cast_integrator_rules = []
    pylink_qmc_instantiate_make_amplitudes_rules = []
    pylink_qmc_instantiate_make_integral_rules = []
    pylink_qmc_instantiate_amplitude_integral_rules = []
    pylink_qmc_transform_instantiation_rules = []
    pylink_qmc_extern_rules = []
    pylink_qmc_case_rules = []
    for pylink_qmc_transform in pylink_qmc_transforms:
        pylink_qmc_dynamic_cast_integrator_rules.append(pylink_qmc_dynamic_cast_integrator[pylink_qmc_transform])
        pylink_qmc_instantiate_make_amplitudes_rules.append(pylink_qmc_instantiate_make_amplitudes[pylink_qmc_transform])
        pylink_qmc_instantiate_make_integral_rules.append(pylink_qmc_instantiate_make_integral[pylink_qmc_transform])
        pylink_qmc_instantiate_amplitude_integral_rules.append(pylink_qmc_instantiate_amplitude_integral[pylink_qmc_transform])
        pylink_qmc_extern_rules.append(pylink_qmc_extern_translation[pylink_qmc_transform])
        pylink_qmc_case_rules.append(pylink_qmc_case_translation[pylink_qmc_transform])
    for i, pylink_qmc_transform in enumerate(chunks(pylink_qmc_transforms, 5)):
        pylink_qmc_transform_instantiation_rules.append(
            {
                'src': 'qmc_template_instantiations.cpp',
                'dest': 'qmc_template_instantiations_%i.cpp' % i,
                'replacements': {
                    'pylink_qmc_transforms': ' '.join([pylink_qmc_instantiations_translation[x] for x in pylink_qmc_transform])
                }
            }
        )

    # check if complex return type required
    enforce_complex = any(g.enforce_complex for g in package_generators)
    contour_deformation = any(g.contour_deformation_polynomial for g in package_generators)
    need_complex = len(complex_parameters) > 0 or enforce_complex or contour_deformation

    replacements_in_files = {
                                'name' : name,
                                'name_as_char_pack' : ','.join([("'" + char + "'") for char in name]),
                                'number_of_integrals' : len(package_generators),
                                'number_of_real_parameters' : len(real_parameters),
                                'names_of_real_parameters' : make_cpp_list(real_parameters),
                                'number_of_complex_parameters' : len(complex_parameters),
                                'names_of_complex_parameters' : make_cpp_list(complex_parameters),
                                'have_complex_parameters': len(complex_parameters) > 0,
                                'number_of_amplitudes' : len(sums),
                                'integral_names' : sub_integral_names,
                                'integral_initialization' : integral_initialization,
                                'integral_initialization_with_contour_deformation': integral_initialization_with_contour_deformation,
                                'weighted_integral_includes' : weighted_integral_includes,
                                'weighted_integral_sum_initialization' : weighted_integral_sum_initialization,
                                'number_of_regulators' : len(regulators),
                                'names_of_regulators' : make_cpp_list(regulators),
                                'requested_orders' : ','.join(map(str,requested_orders)),
                                'pySecDec_version' : version,
                                'python_executable' : sys.executable,
                                'python_version' : sys.version,
                                'pySecDec_git_id' : git_id,
                                'contrib_dirname' : pySecDecContrib.dirname,
                                'date_time' : strftime("%a %d %b %Y %H:%M"),
                                'nested_series_type' : nested_series_type,
                                'pylink_qmc_dynamic_cast_integrator' : '\n        '.join(pylink_qmc_dynamic_cast_integrator_rules),
                                'pylink_qmc_instantiate_make_amplitudes': '\n    '.join(pylink_qmc_instantiate_make_amplitudes_rules),
                                'pylink_qmc_instantiate_make_integral': '\n        '.join(pylink_qmc_instantiate_make_integral_rules),
                                'pylink_qmc_instantiate_amplitude_integral': '\n        '.join(pylink_qmc_instantiate_amplitude_integral_rules),
                                'pylink_qmc_externs': ' '.join(pylink_qmc_extern_rules),
                                'pylink_qmc_cases': ' '.join(pylink_qmc_case_rules),
                                'need_complex': int(bool(need_complex)),
                                'enforce_complex': int(bool(enforce_complex)),
                                'enforce_complex_return_type': int(bool(enforce_complex)),  # make sure that this is either ``0`` or ``1``
                                'contour_deformation': int(bool(contour_deformation))
    }
    filesystem_replacements = {
                                  'integrate_name.cpp' : 'integrate_' + name + '.cpp',

                                  # the following files are integral specific
                                  'name_weighted_integral.cpp' : None,
                                  'name_weighted_integral.hpp' : None,

                                  # the following files are parsed separately
                                  'name.hpp': None,
                                  'qmc_template_instantiations.cpp': None
    }

    # get path to the directory with the template files (path relative to directory with this file: "./templates/")
    template_sources = os.path.join(os.path.split(os.path.abspath(__file__))[0],'templates','sum_package')
    parse_template_tree(template_sources, name, replacements_in_files, filesystem_replacements)

    # Parse pylink files
    for pylink_qmc_transform in pylink_qmc_transform_instantiation_rules:
        pylink_qmc_transform_replacements = replacements_in_files.copy()
        pylink_qmc_transform_replacements.update(pylink_qmc_transform['replacements'])
        parse_template_file(os.path.join(template_sources, 'pylink', pylink_qmc_transform['src']),
                            os.path.join(name,             'pylink', pylink_qmc_transform['dest']),
                            pylink_qmc_transform_replacements)

    original_working_directory = os.getcwd()

    os.mkdir(os.path.join(name, "disteval"))
    os.mkdir(os.path.join(name, "disteval", "coefficients"))
    try:
        os.chdir(name)
        with open("integral_names.txt","w") as f:
            f.write(name+"\n\n"+"\n".join(sub_integral_names.split()))

        # call package generator for every integral
        if processes > 1:
            with Pool(processes) as pool:
                template_replacements = pool.starmap(_generate_one_term, [(
                        j, sums, complex_parameters,
                        name, package_generator._replace(processes=1),
                        pylink_qmc_transforms, real_parameters, regulators,
                        replacements_in_files, requested_orders,
                        template_sources
                    )
                    for j, package_generator in enumerate(package_generators)
                ])
        else:
            template_replacements = [_generate_one_term(
                    j, sums, complex_parameters,
                    name, package_generator,
                    pylink_qmc_transforms, real_parameters, regulators,
                    replacements_in_files, requested_orders,
                    template_sources
                )
                for j, package_generator in enumerate(package_generators)
            ]

        replacements_in_files['number_of_integration_variables'] = max(
            t['number_of_integration_variables']
            for coeffs_lo, coeff_ho, t in template_replacements
        )

    finally:
        os.chdir(original_working_directory)

    # Parse sum_package header file
    parse_template_file(os.path.join(template_sources, 'name.hpp'),  # source
                        os.path.join(name, name + '.hpp'),  # dest
                        replacements_in_files)

    with open(os.path.join(name, "disteval", name + ".json"), "w") as f:
        json.dump({
            "name": name,
            "type": "sum",
            "max_dimension": replacements_in_files["number_of_integration_variables"],
            "regulators": [str(p) for p in regulators],
            "realp": [str(p) for p in real_parameters],
            "complexp": [str(p) for p in complex_parameters],
            "complex_result": bool(need_complex),
            "requested_orders": requested_orders,
            "integrals": [p.name for p in package_generators],
            "sums": {
                sum_name : [
                    {
                        "integral": package_generators[genidx].name,
                        "coefficient": f"{package_generators[genidx].name}_coefficient{sumidx}.txt",
                        "coefficient_lowest_orders": coeffs_lo[sumidx],
                        "coefficient_highest_orders": list(map(int, coeff_ho))
                    }
                    for genidx, (coeffs_lo, coeff_ho, t) in enumerate(template_replacements)
                    if sumidx in coeffs_lo
                ]
                for sumidx, sum_name in enumerate(sums.keys())
            }
        }, f, indent=2)
    # Return template replacements of last integral processed (for 1 integral case this emulates what code_writer.make_package does)
    return template_replacements[-1][2]
