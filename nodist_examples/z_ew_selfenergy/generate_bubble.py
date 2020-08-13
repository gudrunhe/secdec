#! /usr/bin/env python
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    #
    # Example provided by Krzysztof Grzanka 23/01/2020
    #

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ['k1 ** 2',' (k1 - k2) ** 2',' k2 ** 2',' (k1 - k3) ** 2',' (k2 - k3) ** 2',' k3**2 - 1',' (k1 + p1) ** 2',' (k2 + p1) ** 2'],
    loop_momenta = [ 'k1', 'k2', 'k3' ],
    powerlist = [1, 1, 1, 1, 1, 1, 1, 1],
    replacement_rules = [
                            ('p1*p1','1')
                        ]
    )

    Mandelstam_symbols = []
    mass_symbols = []


    loop_package(

    name = 'bubble',
    additional_prefactor = 'exp(3*EulerGamma*eps)',
    loop_integral = li,
    real_parameters = Mandelstam_symbols+mass_symbols,
    complex_parameters = [],
    requested_order = 0,
    form_optimization_level = 2,
    form_work_space = '2G',
    decomposition_method = 'iterative',
    normaliz_executable='normaliz',
    contour_deformation = True,
    split = True,
    )
