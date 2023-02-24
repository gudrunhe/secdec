#!/usr/bin/env python3
import sympy as sp
import json
from pySecDec.integral_interface import DistevalLibrary

if __name__ == "__main__":

    with open('muon_decay2L_mass_expansion/prefactor_exponent_by_name.json') as f:
        prefactor_exponent_by_name = json.load(f)

    final_result = ''
    #integrate and format outputs
    for name, exponent in prefactor_exponent_by_name.items():
        loop_integral = DistevalLibrary('muon_decay2L_mass_expansion/{0}/disteval/{0}.json'.format(name))
        str_result = loop_integral(parameters={'s' : 3, 'mwsq' : 0.78, 'mzsq' : 1.0}, verbose=True)
        result = sp.sympify(str_result)
        value = result[0].subs({"plusminus": 0})
        error = result[0].coeff("plusminus")

        mass_factor = 'mtsq^(' + exponent + ')'
        eps0 = '(' + str(value.coeff('eps',0)) + ') +/- (' + str(error.coeff('eps',0)) + ')'
        eps1 = '(' + str(value.coeff('eps',-1)) + ') +/- (' + str(error.coeff('eps',-1)) + ')'
        eps2 = '(' + str(value.coeff('eps',-2)) + ') +/- (' + str(error.coeff('eps',-2)) + ')'

        final_result += mass_factor + '*{\n' + '   eps^0[' + eps0 + '] +' + '\n' + '   eps^(-1)[' + eps1 + '] +' + '\n' + '   eps^(-2)[' + eps2 + ']}' + ' +' + '\n'

    final_result = final_result[:len(final_result)-2] #Remove final plus sign and linebreak

    #print result
    print('Result:','\n',final_result)