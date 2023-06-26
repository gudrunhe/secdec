#! /usr/bin/env python
from pySecDec.code_writer import make_package
import pySecDec as psd

if __name__ == "__main__":

    #
    # Original Integral
    #
    #li = psd.loop_integral.LoopIntegralFromPropagators(
    #propagators = ['k1^2', '(k1+p1+p2)^2', '(k1-k2)^2', '(k1-k2+p1)^2-mZ^2', '(k2)^2', '(k2+p2)^2'],
    #loop_momenta = ['k1','k2'],
    #
    #external_momenta = ['p1','p2','p3'],
    #
    #replacement_rules = [
    #                        ('p1*p1', 0),
    #                        ('p2*p2', 0),
    #                        ('p3*p3', 's'),
    #                        ('p1*p2', 's/2'),
    #                        ('p2*p3', '-s/2'),
    #                        ('p1*p3', '-s/2'),
    #                        ('s', 'mZ^2'),
    #                        ('mZ', 1)
    #                    ]
    #
    #)
    #
    #Mandelstam_symbols = []
    #mass_symbols = []
    #
    #loop_package(
    #name = 'triangle2L_wu',
    #additional_prefactor = '-exp(2*EulerGamma*eps)',
    #loop_integral = li,
    #real_parameters = Mandelstam_symbols + mass_symbols,
    ## the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    #requested_order = 0,
    ## the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    #form_optimization_level = 2,
    ## the WorkSpace parameter for FORM
    #form_work_space = '1G',
    ## the method to be used for the sector decomposition
    ## valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
    #decomposition_method = 'iterative',
    ## if you choose ``geometric[_ku]`` and 'normaliz' is not in your
    ## $PATH, you can set the path to the 'normaliz' command-line
    ## executable here
    ##normaliz_executable='/path/to/normaliz',
    ## whether or not to produce code to perform the contour deformation
    ## contour deformation is not required if we only want to compute euclidean points (all Mandelstam invariants negative)
    #contour_deformation = True,
    ## there are singularities at one due to ``s = mZ^2``
    #split = True,
    #)


    #
    # Rescaled Integral
    #
    Mandelstam_symbols = []
    mass_symbols = []

    feynman_parameters = ['x0', 'x1', 'x2', 'x3', 'x4', 'x5']
    
    # Note: we could also introduce rescale factors r0, r1 etc which are real_parameters, and set them during integration (might generate larger code)
    rescale_factors = ['(63/73)', '(29/31)', '(27/37)', '(53/59)', '(81/97)', '(20/23)'] # rescale all variables by ratios of 2 digit primes
    #rescale_factors = ['(63/73)', '(29/31)', '(1)', '(53/59)', '(1)', '(1)'] # just rescale 3 variables
    jacobian_factor = '*'.join(rescale_factors)
    replacements = [ (parameter, '('+factor+'*'+parameter+')' ) for parameter,factor in zip(feynman_parameters,rescale_factors) ]

    # li.Gamma_factor
    # li.exponentiated_U
    # li.exponentiated_F
    additional_prefactor = '-exp(2*EulerGamma*eps)'
    gamma_factor ='gamma(2*eps + 2)'
    exponentiated_U = '( + (1)*x3*x5 + (1)*x3*x4 + (1)*x2*x5 + (1)*x2*x4 + (1)*x1*x5 + (1)*x1*x4 + (1)*x1*x3 + (1)*x1*x2 + (1)*x0*x5 + (1)*x0*x4 + (1)*x0*x3 + (1)*x0*x2)**(3*eps)'
    exponentiated_F = '( + (1)*x3**2*x5 + (1)*x3**2*x4 + (1)*x2*x3*x5 + (1)*x2*x3*x4 + (1)*x1*x3*x5 + (1)*x1*x3*x4 + (1)*x1*x3**2 + (-1)*x1*x2*x4 + (1)*x1*x2*x3 + (1)*x0*x3*x4 + (1)*x0*x3**2 + (1)*x0*x2*x3 + (-1)*x0*x1*x5 + (-1)*x0*x1*x4 + (-1)*x0*x1*x3 + (-1)*x0*x1*x2)**(-2*eps - 2)'

    # Note: DANGER!!! would replace x11 -> ((1)*x1)1
    for parameter,replacement in replacements:
        exponentiated_U = exponentiated_U.replace(parameter,replacement)
        exponentiated_F = exponentiated_F.replace(parameter,replacement)

    print(exponentiated_U)
    print(exponentiated_F)

    make_package(

    name = 'triangle2L_wu',

    integration_variables = ['x%i' % i for i in range(len(feynman_parameters))],
    regulators = ['eps'],

    requested_orders = [0],

    polynomials_to_decompose = [exponentiated_U,exponentiated_F],
    polynomial_names = ['U','F'],
    other_polynomials = [jacobian_factor],

    prefactor = '('+additional_prefactor+')'+'*'+'('+gamma_factor+')',
    real_parameters = Mandelstam_symbols+mass_symbols,

#    form_optimization_level = 2,
#    form_work_space = '500M',

    contour_deformation_polynomial = 'F',
    positive_polynomials = 'U',
    decomposition_method = 'iterative',

    )
