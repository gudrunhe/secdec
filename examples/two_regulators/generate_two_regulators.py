#! /usr/bin/env python
from pySecDec import make_package

make_package(

name='two_regulators',
integration_variables = ['z%i' % i for i in range(2)],

# the order here defines the order of the expansion
regulators = ['alpha', 'eps'],

#complex_parameters = ['A'],
#functions = ['F1', 'F2'],

prefactor = '''exp(-EulerGamma*(2*eps+alpha))''',

polynomials_to_decompose = ['z0**(-1-2*eps-alpha)','(1-z0)**(-1+2*eps+alpha)','z1**(-1+alpha/2)'],
remainder_expression = 'exp(-z0/(1-z0))',

# the highest orders of the final regulator expansion
# the order here matches the order of ``regulators``
requested_orders = [0, 0],

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '500M',


)

# analytic result = 2/alpha*gamma(-2*eps-alpha)*exp(-EulerGamma*(2*eps+alpha));
