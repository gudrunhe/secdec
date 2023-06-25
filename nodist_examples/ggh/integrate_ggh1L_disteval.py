#!/usr/bin/env python3
from pySecDec.integral_interface import DistevalLibrary, IntegralLibrary
import sympy as sp

if __name__ == "__main__":

    # load the library
    dlib = DistevalLibrary('ggh1L/disteval/ggh1L.json')

    # integrate
    #dlib_str_result1 = dlib(parameters={"s12": 4.523321402575598, 'mt2': 1.0}, verbose=True, epsrel = 1e-7, epsabs = 1e-10, lattice_candidates=0, format='mathematica')
    dlib_str_result1 = dlib(parameters={"s12": 2.11426927144226, 'mt2': 1.0}, verbose=True, epsrel = 1e-7, epsabs = 1e-10, lattice_candidates=0, format='mathematica')

    print(dlib_str_result1)

