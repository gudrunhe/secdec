#! /usr/bin/env python
import shutil
from pySecDec import make_package

name = 'userdefined_cpp'

make_package(

name=name,
integration_variables = ['z%i' % i for i in range(4)],

# the order here defines the order of the expansion
regulators = ['eps'],
real_parameters = ['alpha'],
functions = ['dum1', 'dum2'],
#
polynomials_to_decompose = ['(z0+z1)**(-2-2*eps)', 'z2**(-1-4*eps)'],
remainder_expression = '(dum1(z0,z1,z2,z3) + 5*eps*z0)**(1+eps) * dum2(z0,z1,alpha)**(2-6*eps)',
#remainder_expression = '(z0**2+z1**3+z2**4+z3**5 + 4*z0*z1*z2*z3+2-z0**2*z1**3*z2**4*z3**5 + 5*eps*z0)**(1+eps) * (z0**2 + z1**2 +alpha**2 + 4*z0*z1+3*z0**2*z1**2 - sqrt(z0*z1*alpha))**(2-6*eps)',


# the highest orders of the final regulator expansion
# the order here matches the order of ``regulators``
requested_orders = [1],

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '500M',

)

# check the generated "functions.hpp"
with open(name + '/src/functions.hpp') as generated_header_file:
    generated_header = generated_header_file.read()
with open('functions_template_ordering_1.hpp') as target_header_1:
    with open('functions_template_ordering_2.hpp') as target_header_2:
        assert generated_header == target_header_1.read() or generated_header == target_header_2.read(), 'mismatch between generated and expected "functions.hpp"'

# copy 'functions.hpp' (predefined for this example) to required directory
shutil.copy('functions_implementation.hpp',name+'/src/functions.hpp')
