#! /usr/bin/env python
import shutil
from pySecDec import make_package

make_package(

name='theta_function',
integration_variables = ['z%i' % i for i in range(3)],

# the order here defines the order of the expansion
regulators = ['eps'],
real_parameters = ['delt'],
functions = ['cut1'],

polynomials_to_decompose = ['(z0+z1)**(-2-2*eps)', 'z2**(-1-4*eps)'],
remainder_expression = 'cut1(z1,delt)',


# the highest orders of the final regulator expansion
# the order here matches the order of ``regulators``
requested_orders = [1],

# the optimization level to use in FORM (can be 0, 1, 2, 3)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '500M',


)

# copy 'functions.hpp' (predefined for this example) to required directory
shutil.copy('functions_theta_sav.hpp','theta_function/src/functions.hpp')
