#! /usr/bin/env python
from pySecDec import make_package

# example is integral representation of Hypergeometric function 5F4

make_package(

name='hypergeo5F4',
integration_variables = ['z%i' % i for i in range(4)],

# the order here defines the order of the expansion
regulators = ['eps'],
real_parameters = ['b'],
#functions = [],
prefactor = '''
                    gamma(2*eps)*gamma(4*eps)*gamma(6*eps)*gamma(8*eps)/
                    (gamma(-eps)*gamma(-3*eps)*gamma(-5*eps)*gamma(-7*eps))/
                    gamma(3*eps)/gamma(7*eps)/gamma(11*eps)/gamma(15*eps)
            ''',


#prefactor = '947.4609375*eps**4 - 286765.9822507143*eps**6 + 1.4350164700988792*^6*eps**7 + 2.668449372562483*^7*eps**8',

polynomials_to_decompose = ['z0**(-1-7*eps)','z1**(-1-5*eps)','z2**(-1-3*eps)','z3**(-1-eps)', '(1-z0)**(-1+15*eps)','(1-z1)**(-1+11*eps)','(1-z2)**(-1+7*eps)','(1-z3)**(-1+3*eps)','(1-b*z0*z1*z2*z3)**(-eps)'],
#remainder_expression = '',

# the highest orders of the final regulator expansion
# the order here matches the order of ``regulators``
requested_orders = [4],

# the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
form_optimization_level = 2,

# the WorkSpace parameter for FORM
form_work_space = '500M',

split = True,

)

# analytic result for ``b=0.5`` (obtained using HypExp, T.Huber, D.Maitre, hep-ph/0507094)

#result_analytic = 1 + 0.18953243218436003*eps - 2.2990427423820186*eps^2 + 55.4690190360555*eps^3 - 1014.392422652346*eps^4 + O(eps^5);
