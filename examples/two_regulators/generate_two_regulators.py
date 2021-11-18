#!/usr/bin/env python3
from pySecDec import make_package

if __name__ == "__main__":

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
    
    # Note:
    # * this integral is numerically regulated for z0->1 (see remainder_expression)
    # * default korobov3x3 transforms the integrand strongly close to z0->1
    # * applying korobov3x3 makes the integrand too numerically unstable to integrate (often get nan)
    # We instead use the asymmetric korobov3x1 to avoid introducing too much noise for z0->1
    pylink_qmc_transforms=['korobov3x1']

    )

    # analytic result = 2/alpha*gamma(-2*eps-alpha)*exp(-EulerGamma*(2*eps+alpha));
