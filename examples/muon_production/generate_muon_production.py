import pySecDec as psd
from pySecDec import LoopPackage
from pySecDec import sum_package
import numpy as np

if __name__ == "__main__":
    
    coeffs = [
        '8*((4 - 2*eps)^2*s^2 - 3*(4 - 2*eps)*t^2 - 10*(4 - 2*eps)*t*u - 3*(4 - 2*eps)*u^2 + 2*t^2 + 8*t*u + 2*u^2)/((3 - 2*eps)*s^2)',
        '2*((4 - 2*eps)^3*s^2 - 11*(4 - 2*eps)^2*t^2 + 26*(4 - 2*eps)*t^2 - 24*(4 - 2*eps)^2*t*u + 78*(4 - 2*eps)*t*u - 9*(4 - 2*eps)^2*u^2 + 20*(4 - 2*eps)*u^2 - 16*t^2 - 56*t*u - 12*u^2)/((3 - 2*eps)*s^2)',
        '((4 - 2*eps)*t + 3*(4 - 2*eps)*u - 4*t - 8*u)/s',
        '-(3*(4 - 2*eps)*t + (4 - 2*eps)*u - 8*t - 4*u)/s',
        '-2*((4 - 2*eps)^2*s^2 - 5*(4 - 2*eps)*s^2 + 4*(4 - 2*eps)*u^2 + 4*t^2 + 8*t*u)/((3 - 2*eps)*s)',
        '-t*((4 - 2*eps)^2*s + 9*(4 - 2*eps)*t + (4 - 2*eps)*u - 8*t)/(2*(3 - 2*eps)*s)',
        'u*((4 - 2*eps)^2*s + (4 - 2*eps)*t + 9*(4 - 2*eps)*u - 8*u)/(2*(3 - 2*eps)*s)',
        '-2*((4 - 2*eps)^2*s^2 - 5*(4 - 2*eps)*s^2 + 4*(4 - 2*eps)*u^2 + 4*t^2 + 8*t*u)/((3 - 2*eps)*s)',
        '-t*((4 - 2*eps)^2*s + 9*(4 - 2*eps)*t + (4 - 2*eps)*u - 8*t)/(2*(3 - 2*eps)*s)',
        'u*((4 - 2*eps)^2*s + (4 - 2*eps)*t + 9*(4 - 2*eps)*u - 8*u)/(2*(3 - 2*eps)*s)',
        '-t*(3*(4 - 2*eps)^2*s^2 - 3*(4 - 2*eps)*t^2 - 30*(4 - 2*eps)*t*u - 11*(4 - 2*eps)*u^2 + 24*t*u + 8*u^2)/(2*(3 - 2*eps)*s)',
        'u*(3*(4 - 2*eps)^2*s^2 - 11*(4 - 2*eps)*t^2 - 30*(4 - 2*eps)*t*u - 3*(4 - 2*eps)*u^2 + 8*t^2 + 24*t*u)/(2*(3 - 2*eps)*s)']

    N_coeffs = [
        '0',
        '2*(-(4 - 2*eps)^2*s^2 + 4*(4 - 2*eps)*t^2 + 12*(4 - 2*eps)*t*u + 4*(4 - 2*eps)*u^2 - 4*t^2 - 16*t*u - 4*u^2)/((3 - 2*eps)*s^2)',
        '0',
        '0',
        '0',
        '0',
        '0',
        '0',
        '0',
        '0',
        '0',
        '0']

    additional_prefactor = 'gamma(1-2*eps)/(gamma(1-eps)*gamma(1-eps)*gamma(1 + eps))'
    
    def B0(p_sq, name):
        li = psd.LoopIntegralFromGraph(
                internal_lines = [[0,[1,2]],[0,[2,1]]],
                external_lines = [['p',1],['p',2]],
                replacement_rules = [('p*p', p_sq)])
        real_parameters = []
        if not p_sq == 0: #Pass the momentum as a symbolic parameter if it is not 0
            real_parameters.append(p_sq)
        return LoopPackage(name, loop_integral = li, real_parameters = real_parameters, 
                            decomposition_method = 'geometric', requested_orders = [0], additional_prefactor = additional_prefactor)

    def C0(p1_sq, p2_sq, p12_sq, name):
        li = psd.LoopIntegralFromGraph(
                internal_lines = [[0,[1,2]],[0,[2,3]],[0,[3,1]]],
                external_lines = [['p1',1],['p2',2],['p3',3]],
                replacement_rules = [
                                    ('p1*p1', p1_sq),
                                    ('p2*p2', p2_sq),
                                    ('p3*p3', p12_sq)
                                    ])
        real_parameters = [p for p in [p1_sq, p2_sq, p12_sq] if p != 0] #Pass the momenta as symbolic parameters if they are not 0
        return LoopPackage(name, loop_integral = li, real_parameters = real_parameters, 
                            decomposition_method = 'geometric', requested_orders = [0], additional_prefactor = additional_prefactor)

    def D0(p1_sq, p2_sq, p3_sq, p4_sq, p12_sq, p23_sq, name):
        li = psd.LoopIntegralFromGraph(
                internal_lines = [[0,[1,2]],[0,[2,3]],[0,[3,4]],[0,[4,1]]],
                external_lines = [['p1',1],['p2',2],['p3',3], ['p4',4]],
                replacement_rules = [('p1*p1', p1_sq),
                                                ('p2*p2', p2_sq),
                                                ('p3*p3', p3_sq),
                                                ('p4*p4', p4_sq),
                                                ('p3*p2', str(p23_sq) + '/2' + '-' + str(p2_sq) + '/2' + '-' + str(p3_sq) + '/2'),
                                                ('p1*p2', str(p12_sq) + '/2' + '-' + str(p1_sq) + '/2' + '-' + str(p2_sq) + '/2'),
                                                ('p1*p4', str(p23_sq) + '/2' + '-' + str(p1_sq) + '/2' + '-' + str(p4_sq) + '/2'),
                                                ('p2*p4', '-' + str(p12_sq) + '/2' + '-' + str(p23_sq) + '/2' + '-' + str(p2_sq) + '/2' + '-' + str(p4_sq) + '/2'),
                                                ('p1*p3', '-' + str(p12_sq) + '/2' + '-' + str(p23_sq) + '/2' + '-' + str(p1_sq) + '/2' + '-' + str(p3_sq) + '/2'),
                                                ('p3*p4', str(p12_sq) + '/2' + '-' + str(p3_sq) + '/2' + '-' + str(p4_sq) + '/2')
                                                ])
        real_parameters = [p for p in [p1_sq, p2_sq, p3_sq, p4_sq, p12_sq, p23_sq] if p != 0] #Pass the momenta as symbolic parameters if they are not 0
        return LoopPackage(name, loop_integral = li, real_parameters = real_parameters, 
                            decomposition_method = 'geometric', requested_orders = [0], additional_prefactor = additional_prefactor)

    integrals = [
        B0(0, 'B00'),
        B0('s', 'B0s'),
        B0('t', 'B0t'),
        B0('u', 'B0u'),
        C0(0, 0, 's', 'C00s'),
        C0(0, 0, 't', 'C00t'),
        C0(0, 0, 'u', 'C00u'),
        C0(0, 's', 0, 'C0s0'),
        C0(0, 't', 0, 'C0t0'),
        C0(0, 'u', 0, 'C0u0'),
        D0(0,0,0,0, 's', 't', 'D0000st'),
        D0(0,0,0,0, 's', 'u', 'D0000su')
    ]

    sum_package(
            'full_amplitude',
            integrals,
            coefficients = {'O(1)': coeffs, 'O(Nf)': N_coeffs},
            regulators = ['eps'],
            requested_orders = [0],
            real_parameters = ['s', 't', 'u'], #Make sure that the list of real parameters contain every symbolic kinematic invariant defined in 'all_integrals'
            complex_parameters = [],
            processes = 30)
