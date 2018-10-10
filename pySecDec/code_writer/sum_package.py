from ..metadata import version, git_id
from .template_parser import parse_template_file, parse_template_tree
from ..misc import sympify_symbols, make_cpp_list
from ..algebra import Polynomial, ExponentiatedPolynomial
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

            Symbol SecDecInternalsDUMMY;
            Symbols `regulators';
            Symbols `parameters';

            %(expression_definitions)s
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
