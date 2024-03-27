#!/usr/bin/env python3
import shutil
from itertools import permutations
from pySecDec import make_package

if __name__ == "__main__":

    name = 'userdefined_cpp'

    make_package(

    name=name,
    integration_variables = ['x','y','z'],

    # use integration by parts rather than subtraction for ``x`` and ``y``
    # generate subtraction term for ``z`` to avoid derivative of HeavisideTheta
    ibp_power_goal = [0,0,-1],

    # the order here defines the order of the expansion
    regulators = ['eps'],
    functions = ['HeavisideTheta','func'],

    polynomials_to_decompose = ['(x+y)**(-2+eps)','z**(-1+eps)'],
    remainder_expression = 'HeavisideTheta(1/4-z)*func(y)',
    # full analytic result for func(y)=1:
    #   (4^-eps (-2 + 2^eps))/((-1 + eps) eps^2)
    # = 1/eps^2 + (1 - Log[2] - Log[4])/eps +  1/2 (2 - 2 Log[2] - Log[2]^2 - 2 Log[4] + 2 Log[2] Log[4] + Log[4]^2) + O[eps]
    # = 1.0000000000000000000/eps^2 - 1.0794415416798359283/eps + 0.60214400703386905808 + O[eps]

    # the highest orders of the final regulator expansion
    # the order here matches the order of ``regulators``
    requested_orders = [0],

    decomposition_method='iterative_no_primary'
    )

    # check the generated "functions.hpp"
    cpp_function_declarations = (
    '    template<typename T0>\n    integrand_return_t HeavisideTheta(T0 arg0);',
    '    template<typename T0>\n    integrand_return_t dfuncd0(T0 arg0);',
    '    template<typename T0>\n    integrand_return_t func(T0 arg0);'
    )
    with open(name + '/' + name + '_integral/src/functions.hpp') as generated_header_file:
        generated_header = generated_header_file.read()
    with open('functions_template.hpp') as target_header_template_file:
        target_header_template = target_header_template_file.read()
    matched_header = False
    for ordering in permutations(cpp_function_declarations):
        if generated_header == target_header_template % ordering:
            matched_header = True
            break
    assert matched_header, 'mismatch between generated and expected "functions.hpp"'

    # copy 'functions.hpp' (predefined for this example) to required directory
    shutil.copy('functions_implementation.hpp',name + '/' + name + '_integral/src/functions.hpp')
