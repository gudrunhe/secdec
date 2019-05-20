from ..metadata import version, git_id
from .template_parser import parse_template_file, parse_template_tree
from ..misc import sympify_symbols, make_cpp_list
from ..loop_integral.loop_package import loop_package

from time import strftime
from itertools import repeat
import os, sys, shutil, subprocess
import numpy as np
import sympy as sp

class Coefficient(object):
    r'''
    Store a coefficient expressed as a product of terms
    in the numerator and a product of terms in the
    denominator.

    :param numerators:
        iterable of str;
        The terms in the numerator.

    :param denominators:
        iterable of str;
        The terms in the denominator.

    :param regulators:
        iterable of strings or sympy symbols;
        The symbols denoting the regulators.

    :param parameters:
        iterable of strings or sympy symbols;
        The symbols other parameters.

        '''
    def __init__(self, numerators, denominators, regulators, parameters):
        for input_value, input_name in \
            zip(
                    ( numerators ,  denominators ,  regulators ,  parameters ),
                    ('numerators', 'denominators', 'regulators', 'parameters')
               ):
            assert np.iterable(input_value), "`%s` must be iterable." % input_name
            if input_name in ('numerators', 'denominators'):
                for item in input_value:
                    assert isinstance(item, str), "All `%s` must be strings." % input_name

        self.numerators = list(numerators)
        self.denominators = list(denominators)
        self.regulators = sympify_symbols(regulators, 'All `regulators` must be symbols.')
        self.parameters = sympify_symbols(parameters, 'All `parameters` must be symbols.')

    def process(self, form=None, workdir='form_tmp', keep_workdir=False):
        r'''
        Calculate and return the lowest orders of the coefficients in
        the regulators and a string defining the expressions "numerator",
        "denominator", and "regulator_factor".

        :param form:
            string or None;
            If given a string, interpret that string as the command
            line executable `form`. If ``None``, try
            ``$SECDEC_CONTRIB/bin/form`` and, if the
            environment variable ``$SECDEC_CONTRIB`` is not set,
            ``form``.

        :param workdir:
            string;
            The directory for the communication with `form`.
            A directory with the specified name will be created
            in the current working directory. If the specified
            directory name already exists, an :class:`OSError`
            is raised.

            .. note::
                The communication with `form` is done via
                files.

        :param keep_workdir:
            bool;
            Whether or not to delete the `workdir` after execution.

        '''
        # get form command-line executable
        if form is None:
            try:
                # use "$SECDEC_CONTRIB/bin/form" if "$SECDEC_CONTRIB" is defined
                form = os.path.join(os.environ['SECDEC_CONTRIB'], 'bin', 'form')
            except KeyError:
                # "$SECDEC_CONTRIB" is not defined --> let the system find "form"
                form = 'form'
        else:
            assert isinstance(form, str), "`form` must be a string."

        # define the form program to run
        program = '''
            #Define OUTFILE "%(outfile)s"
            #Define regulators "%(regulators)s"
            #Define parameters "%(parameters)s"
            #Define numberOfnum "%(number_of_num)i"
            #Define numberOfden "%(number_of_den)i"

            Symbols SecDecInternalsDUMMY,`regulators';
            Symbols I,`parameters';

            %(expression_definitions)s

            * convert I to i_
            multiply replace_(I,i_);
            .sort:read;

            * factor out the global power of each regulator from each numerator and denominator
            #Do regulator = {`regulators',}
              #If x`regulator' != x
                #Do nORd = {num,den}
                  #Do idx = 1,`numberOf`nORd''
                    skip; nskip `nORd'`idx';
                    #$min`regulator'`nORd'`idx' = maxpowerof_(`regulator');
                    if( match(`regulator'^SecDecInternalsDUMMY?$this`regulator'pow) );
                        if($min`regulator'`nORd'`idx' > $this`regulator'pow) $min`regulator'`nORd'`idx' = $this`regulator'pow;
                    else;
                        if($min`regulator'`nORd'`idx' > 0) $min`regulator'`nORd'`idx' = 0;
                    endif;
                    ModuleOption,local,$this`regulator'pow;
                    ModuleOption,minimum,$min`regulator'`nORd'`idx';
                    .sort:min`nORd'`idx';
                    skip; nskip `nORd'`idx';
                    multiply `regulator'^(-(`$min`regulator'`nORd'`idx''));
                    .sort:mul`nORd'`idx';
                  #EndDo
                #EndDo

                #$global`regulator'fac =
                  #Do idx = 1,`numberOfnum'
                    + $min`regulator'num`idx'
                  #EndDo
                  #Do idx = 1,`numberOfden'
                    - $min`regulator'den`idx'
                  #EndDo
                ;

              #EndIf
            #EndDo

            * convert i_ to I
            multiply replace_(i_,I);
            .sort:beforeWrite;

            #write <`OUTFILE'> "numerator = 1"
            #Do idx = 1,`numberOfnum'
              #write <`OUTFILE'> "*(%%E)", num`idx'
            #EndDo
            #write <`OUTFILE'> ";"

            #write <`OUTFILE'> "denominator = 1"
            #Do idx = 1,`numberOfden'
              #write <`OUTFILE'> "*(%%E)", den`idx'
            #EndDo
            #write <`OUTFILE'> ";"

            #write <`OUTFILE'> "regulator_factor = 1"
            #Do regulator = {`regulators',}
              #If x`regulator' != x
                #write <`OUTFILE'> "*`regulator'^(`$global`regulator'fac')"
              #EndIf
            #EndDo
            #write <`OUTFILE'> ";"

            .end

        '''.replace('\n            ','\n')
        expression_definitions = []
        for idx, numerator in enumerate(self.numerators):
            expression_definitions.append('Local num%i = %s;' % (idx+1,numerator))
        for idx, denominator in enumerate(self.denominators):
            expression_definitions.append('Local den%i = %s;' % (idx+1,denominator))
        outfile = 'results.txt'
        program = program % dict(
                                     expression_definitions='\n'.join(expression_definitions),
                                     outfile=outfile,
                                     number_of_num=len(self.numerators),
                                     number_of_den=len(self.denominators),
                                     regulators=','.join(map(str,self.regulators)),
                                     parameters=','.join(map(str,self.parameters))
                                )

        # run form program
        os.mkdir(workdir)
        try:
            # dump program to file
            with open(os.path.join(workdir, 'program.frm'), 'w') as f:
                f.write(program)

            # redirect stdout
            with open(os.path.join(workdir, 'stdout'), 'w') as stdout:
                # redirect stderr
                with open(os.path.join(workdir, 'stderr'), 'w') as stderr:
                    # run form
                    #    subprocess.check_call --> run form, block until it finishes and raise error on nonzero exit status
                    try:
                        subprocess.check_call([form, 'program'], stdout=stdout, stderr=stderr, cwd=workdir)
                    except OSError as error:
                        if form not in str(error):
                            error.filename = form
                        raise

            # collect output
            with open(os.path.join(workdir, outfile), 'r') as f:
                form_output = f.read()
            assert len(form_output) > 0, "form output file (" + outfile + ") empty"

        finally:
            if not keep_workdir:
                shutil.rmtree(workdir)

        # read lowest_orders from form output
        lowest_orders = np.empty(len(self.regulators), dtype=int)
        regulator_factor = form_output[form_output.rfind('regulator_factor'):]
        regulator_factor = regulator_factor[:regulator_factor.find(';')].replace('\n','').replace(' ','')
        regulator_factors = regulator_factor.split('=')[1].split('*')[1:]
        assert len(regulator_factors) == len(self.regulators)
        for idx,(regulator,factor) in enumerate(zip(self.regulators,regulator_factors)):
            form_regulator, power = factor.split('^')
            assert form_regulator == str(regulator), '%s != %s' % (form_regulator, regulator)
            lowest_orders[idx] = int(power.lstrip('(').rstrip(')'))

        return lowest_orders, form_output

# TODO: high-level test
def sum_package(name, package_generators, generators_args, regulators, requested_orders,
                real_parameters=[], complex_parameters=[], coefficients=None,
                form_executable=None):
    r'''
    Decompose, subtract and expand an expression of the form

    .. math::
        \sum_i c_i \ \int f_i

    Generate a c++ package with an optmized algorithm to
    evaluate the integrals numerically.

    :param name:
        string;
        The name of the c++ namepace and the output
        directory.

    :param package_generators:
        iterable of functions;
        The generators functions for the integrals
        :math:`\int f_i`

        .. seealso::
            :func:`pySecDec.code_writer.make_package` and
            :func:`pySecDec.loop_integral.loop_package`

    :param generators_args:
        iterable of dictionaries;
        The keyword arguments that are passed to the
        `package_generators`.

        .. seealso::
            :func:`pySecDec.code_writer.make_package` and
            :func:`pySecDec.loop_integral.loop_package`

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
        iterable of :class:`.Coefficient`, optional;
        The coefficients :math:`c_i` of the integrals.
        :math:`c_i = 1` is assumed if not provided.

    :param form_executable:
        string or None, optional;
        The path to the form exectuable. The argument is passed
        to :meth:`.Coefficient.process`.
        Default: ``None``.

    '''
    print('running "sum_package" for ' + name)

    # convert input iterables to lists
    package_generators = list(package_generators)
    generators_args = list(generators_args)
    regulators = list(regulators)
    requested_orders = list(requested_orders)
    real_parameters = list(real_parameters)
    complex_parameters = list(complex_parameters)

    # prepare coefficients
    empty_coefficient = Coefficient((), (), regulators, ())
    if coefficients is None:
        coefficients = repeat(empty_coefficient)
    else:
        coefficients = [empty_coefficient if c is None else c for c in coefficients]
        assert len(coefficients) == len(package_generators), \
            "`coefficients` must either be ``None`` have the same length as `package_generators`."

    assert len(package_generators) == len(generators_args), \
        "`package_generators` (%i) and `generators_args` (%i) must have the same length." % \
        (len(package_generators), len(generators_args))

    # construct the c++ type "nested_series_t"
    # for two regulators, the resulting code should read:
    # "secdecutil::Series<secdecutil::Series<T>>"
    nested_series_type = 'secdecutil::Series<' * len(requested_orders) + 'T' + '>' * len(requested_orders)

    # define required listings of contributing integrals
    sub_integral_names = []
    weighted_integral_includes = []
    weighted_integral_sum_initialization = []
    for package_generator,generator_args in zip(package_generators,generators_args):
        sub_name = generator_args["name"]
        sub_integral_names.append(sub_name)
        weighted_integral_includes.append( '#include "' + sub_name + '_weighted_integral.hpp"')
        weighted_integral_sum_initialization.append( 'amplitude += %s::make_weighted_integral(real_parameters, complex_parameters);' % sub_name )
    sub_integral_names = ' '.join(sub_integral_names)
    weighted_integral_includes = '\n'.join(weighted_integral_includes)
    weighted_integral_sum_initialization[0] = weighted_integral_sum_initialization[0].replace('+',' ',1)
    weighted_integral_sum_initialization = '\n        '.join(weighted_integral_sum_initialization)

    replacements_in_files = {
                                'name' : name,
                                'number_of_integrals' : len(package_generators),
                                'number_of_real_parameters' : len(real_parameters),
                                'names_of_real_parameters' : make_cpp_list(real_parameters),
                                'number_of_complex_parameters' : len(complex_parameters),
                                'names_of_complex_parameters' : make_cpp_list(complex_parameters),
                                'integral_names' : sub_integral_names,
                                'weighted_integral_includes' : weighted_integral_includes,
                                'weighted_integral_sum_initialization' : weighted_integral_sum_initialization,
                                'number_of_regulators' : len(regulators),
                                'names_of_regulators' : make_cpp_list(regulators),
                                'requested_orders' : ','.join(map(str,requested_orders)),
                                'pySecDec_version' : version,
                                'python_version' : sys.version,
                                'pySecDec_git_id' : git_id,
                                'date_time' : strftime("%a %d %b %Y %H:%M"),
                                'nested_series_type' : nested_series_type
                            }
    filesystem_replacements = {
                                  'integrate_name.cpp' : 'integrate_' + name + '.cpp',
                                  'name.hpp' : name + '.hpp',
                                  'config_name.hpp' : 'config_' + name + '.hpp',

                                  # the following files are integral specific
                                  'name_weighted_integral.cpp' : None,
                                  'name_weighted_integral.hpp' : None
                              }

    # get path to the directory with the template files (path relative to directory with this file: "./templates/")
    template_sources = os.path.join(os.path.split(os.path.abspath(__file__))[0],'templates','sum_package')
    parse_template_tree(template_sources, name, replacements_in_files, filesystem_replacements)

    original_working_directory = os.getcwd()

    try:
        os.chdir(name)

        # call package generator for every integral
        for package_generator, generator_args, coefficient in zip(package_generators, generators_args, coefficients):
            sub_name = generator_args['name']
            replacements_in_files['sub_integral_name'] = sub_name

            generator_args['real_parameters'] = real_parameters
            generator_args['complex_parameters'] = complex_parameters

            # process coefficient
            lowest_coefficient_orders, coefficient_expressions = coefficient.process(form=form_executable)
            with open(os.path.join('lib', sub_name+'_coefficient.txt'),'w') as coeffile:
                coeffile.write(coefficient_expressions)
            replacements_in_files['lowest_coefficient_orders'] = ','.join(map(str,lowest_coefficient_orders))

            if package_generator is loop_package:
                generator_args['requested_order'] = requested_orders[0] - lowest_coefficient_orders[0]
            else:
                generator_args['regulators'] = regulators
                generator_args['requested_orders'] = np.asarray(requested_orders) - lowest_coefficient_orders

            # parse integral specific files
            for ch in 'ch':
                suffix = '_weighted_integral.%spp' % ch
                parse_template_file(os.path.join(template_sources, 'src', 'name' + suffix), # source
                                    os.path.join('src', sub_name + suffix), # dest
                                    replacements_in_files)

            package_generator(**generator_args)

    finally:
        os.chdir(original_working_directory)
