#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(
    loop_momenta = ['k1','k2'],
    Lorentz_indices = ['mu'],
    propagators = ['k1^2 - m1sq', 'k2^2 - m2sq', '(k1+k2)^2 - m3sq'],
    numerator = '2*k1(mu)*k2(mu)',
    )


    loop_package(

    name = 'tadpole2L_rank2',

    loop_integral = li,

    # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
    requested_order = 0,

    real_parameters =['m1sq','m2sq','m3sq'],

    # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
    form_optimization_level = 3,

    # contour deformation is not required because there are no external legs
    contour_deformation = False,

    # make sure the retun type is double
    enforce_complex = True,

    )
