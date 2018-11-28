#! /usr/bin/env python
from pySecDec import make_package

make_package(

name='dummyII',
integration_variables = ['z%i' % i for i in range(4)],

# the order here defines the order of the expansion
regulators = ['eps'],
real_parameters = ['alpha'],
#functions = ['dum1', 'dum2'],
#
polynomials_to_decompose = ['(z0+z1)**(-2-2*eps)', 'z2**(-1-4*eps)'],
other_polynomials = [],
#remainder_expression = '(dum1(z0,z1,z2,z3) + 5*eps*z0)**(1+eps) * dum2(z0,z1,alpha)**(2-6*eps)',
remainder_expression = '(z0**2+z1**3+z2**4+z3**5 + 4*z0*z1*z2*z3+2-z0**2*z1**3*z2**4*z3**5 + 5*eps*z0)**(1+eps) * (z0**2 + z1**2 +alpha**2 + 4*z0*z1+3*z0**2*z1**2 - sqrt(z0*z1*alpha))**(2-6*eps)',


# the highest orders of the final regulator expansion
# the order here matches the order of ``regulators``
requested_orders = [1],

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '500M',


)
