#!/usr/bin/env python3

import pySecDec as psd

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def expand_integral(integral,smallness_parameters):
    common_args = {
        'integration_variables' : integral.integration_variables,
        'regulators' : integral.regulators,
        'requested_orders' : [],
        'expansion_by_regions_order' : 0,
        'polytope_from_sum_of' : [0,1],
        'polynomial_names' : ["U","F"]
    }
    print(bcolors.OKBLUE)
    print('-- ' + name + ' --')
    print('  U: ',integral.exponentiated_U)
    print('  F:',integral.exponentiated_F)
    print('  expanding in:',smallness_parameters[0])
    #print('  measure',integral.measure)
    #print('  numerator',integral.numerator)
    generators_args = psd.make_regions(
        name = name,
        smallness_parameter = smallness_parameters[0],
        polynomials_to_decompose = [integral.exponentiated_U,integral.exponentiated_F],
        **common_args
    )
    print(bcolors.ENDC)
    expand_integrals(generators_args,smallness_parameters,1,common_args=common_args)


def expand_integrals(integrals,smallness_parameters,index,common_args):
    if index == len(smallness_parameters):
        for integral in integrals:
            print(bcolors.OKCYAN)
            print('-- ' + integral.name + ' --')
            print('  U: ',integral.polynomials_to_decompose[0])
            print('  F: ',integral.polynomials_to_decompose[1])
            print('  prefactor: ',integral.prefactor)
            #print('  expanding in: ',smallness_parameters[index])
            print(bcolors.ENDC)
        return
    else:
        for integral in integrals:
            if index == len(smallness_parameters)-1:
                print(bcolors.OKGREEN)
            print('-- ' + integral.name + ' --')
            print('  U: ',integral.polynomials_to_decompose[0])
            print('  F: ',integral.polynomials_to_decompose[1])
            print('  prefactor: ',integral.prefactor)
            print('  expanding in: ',smallness_parameters[index])
            #print(integral)
            try:
                generators_args = psd.make_regions(
                    name = integral.name,
                    smallness_parameter = smallness_parameters[index],
                    polynomials_to_decompose = [integral.polynomials_to_decompose[0],integral.polynomials_to_decompose[1]], # U, F
                     **common_args
                )
            except:
                print('No Expansion')
                continue
            if index == len(smallness_parameters)-1:
                print(bcolors.ENDC)
            expand_integrals(generators_args,smallness_parameters,index+1,common_args=common_args)


if __name__ == "__main__":

    # Formfactor1L
    name = 'formfactor1L_massless_ebr'
    smallness_parameters = ['t']
    internal_lines = [['0',[1,3]],['0',[2,3]],['0',[1,2]]]
    external_lines = [['p1',1],['p2',2],['p3',3]]
    li = psd.LoopIntegralFromGraph(
        internal_lines = internal_lines,
        external_lines = external_lines,
        replacement_rules = [
                            ('p3*p3', '-Qsq'),
                            ('p1*p1', '-t*P1sqr*Qsq'),
                            ('p2*p2', '-t*P2sqr*Qsq'),
                            ('Qsq', 1),
                            ],
    )
    # (Optional) Draw diagrams
    # psd.loop_integral.draw.plot_diagram(internal_lines,external_lines,name)

    expand_integral(li,smallness_parameters)
