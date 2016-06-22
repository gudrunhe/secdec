#! /usr/bin/env python
from pySecDec import make_package

make_package(

target_directory='pySecDec_output',

name='angularities',
integration_variables = ['z%i' % i for i in range(6)],

# the order here defines the order of the expansion
regulators = ['alpha', 'eps'],
complex_parameters = ['A'],
functions = ['DM', 'Z'],
prefactor = '''
                  -1 * 2**6 * exp(EulerGamma)**(-2*eps) / sqrt(pi) * gamma(-4*eps) / gamma(1-eps) / gamma(1/2-eps) * 4 * 2**(-1-2*eps)
                  * 2 * 16 / pi / 2 * (-1) * eps * 2**(-4*eps) * 4**(1+eps) * 64
            ''',
polynomials_to_decompose = ['z3**(-1+4*alpha)', 'z2**(-1-4*eps+2*alpha)', 'z0**(1-4*eps)', '(1-z0) ** (1-4*eps)'], #, '(1+z0) ** (1-4*eps)', '(2-z0**4)**(-1/2-eps)', '(1+z0**2)**(1-4*eps)', '(2-2*z0**4+z0**8)**(-1/2-eps)', '(z5)**(-1-2*eps)', '(2-z5**2)**(-1-eps)', '(z1)**(-1-8*eps)', '(2-z1**4)**(-1-2*eps)', '(2-2*z1**4+z1**8)**(-1-2*eps)', '(z4)**(1-4*eps)', '(1-z1**4*(-2+z1**4)*(2-2*z1**4+z1**8)*z4**4)**(-1/2-eps)', '(1-z1)**(1-4*eps)', '(1+z1)**(1-4*eps)', '(1+z1**2)**(1-4*eps)', '(1+z4**4)**(-3)', '(1+z2**2*(1-z1**4*(2-z1**4)*(2-2*z1**4+z1**8)*(1-z4**4)))**(2*(-1+eps))', '(1+z2**2-z1**4*(2-z1**4)*(2-2*z1**4+z1**8)*(1-z4**4))**(2*eps)'], #, '(z1**2 - 4 * ((1-z1)*z4))**(-1)'],
remainder_expression = 'DM(z1,z2,z3,z4,z5,z0, eps, alpha) * Z(z1,z2,z4, eps)',

# the highest orders of the final regulator expansion
# the order here matches the order of ``regulators``
requested_orders = [0, 0],

# the optimization level to use in FORM (can be 0, 1, 2, 3)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '500M',

# whether or not to increase numerical stability
# Note: This is very extensive concerning both - the algebra and the numerics.
#       It should only be set to ``True`` if numerical instabilities occur.
stabilize = False,



)
