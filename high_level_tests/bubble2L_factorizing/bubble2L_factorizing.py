#!/usr/bin/env python3
from pySecDec.loop_integral import loop_package
import pySecDec as psd

if __name__ == "__main__":

    li = psd.loop_integral.LoopIntegralFromPropagators(
    propagators = ["(k1+p1)^2","(k2-p1)^2","(k2-p2)^2","(k1+p2)^2"],
    loop_momenta = ["k1","k2"],
    replacement_rules = [
                            ("p1^2","0"),
                            ("p2^2","0"),
                            ("p1*p2","es12/2")
                        ],
    regulator = "eps",
    dimensionality = "6-2*eps",
    )

    loop_package(
    name = "bubble2L_factorizing",
    loop_integral = li,
    requested_orders = [-2],
    form_optimization_level = 4,
    decomposition_method = "geometric",
    real_parameters = ['es12'],

    # Contour deformation is not needed but should not alter the result.
    # However, a bug caused wrong results when contour deformation was active.
    contour_deformation = True,

    package_generator=psd.code_writer.make_package
    )
