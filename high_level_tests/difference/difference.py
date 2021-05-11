#!/usr/bin/env python3
from pySecDec import make_package

if __name__ == "__main__":

    make_package(

    name='difference',

    integration_variables = ['z1','z2','z3'],
    regulators = ['eps'],

    polynomials_to_decompose = ['(z1-z2*z3)**(-1+eps)'],

    # the contour deformation adds a "-i * delta" prescription to the polynomial above
    polynomial_names = ['P'],
    contour_deformation_polynomial = 'P',

    # we want to compute up to order "eps**2"
    requested_orders = [2],

    # split the integration region at 1/2 to remap singularities
    # at one to zero
    split = True,

    )


    # analytic result:
    # ( 1.64493406684822643647241516664602518923 + 3.1415926535897932384626433832795028842 * I) * eps ** 0
    # ( 2.08781123053685858754509217178101012328 - 6.2831853071795864769252867665590057684 * I) * eps ** 1
    # (-5.94029019737039970544633397517750766917 + 4.2570651807194096861418776386549427857 * I) * eps ** 2
    # ( 5.77945251635494087034720012662916969501 - 2.2309450542592328953584685107508798034 * I) * eps ** 3 # set ``requested_orders = [3]`` to compute up to order "eps**3"
