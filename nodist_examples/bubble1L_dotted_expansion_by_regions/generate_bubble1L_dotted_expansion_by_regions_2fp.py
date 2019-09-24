# in this example we try to expand the dotted loop using the mass of the particle
# in the loop as the small parameter
# here we do it both ways, expaning in the msq directly, or replacing it by
# z*msq, expanding in z and then using z=1, here we use 2 Feynman parameters
# (i.e. we have two propagators and a powerlist)

import pySecDec as psd
import numpy as np

# for the case where we expand in z
li_m = psd.loop_integral.LoopIntegralFromPropagators(
propagators=("((k+p)**2)", "(k**2-msq_)"),
loop_momenta=["k"],
powerlist=[1,2],
regulators=['eps',],
replacement_rules = [
                        ('p*p', 'psq'),
                        ('msq_', 'msq'),
                    ]
)

li_z = psd.loop_integral.LoopIntegralFromPropagators(
propagators=li_m.propagators, loop_momenta=li_m.loop_momenta, regulators=li_m.regulators, powerlist=li_m.powerlist,
replacement_rules = [
                        ('p*p', 'psq'),
                        ('msq_', 'msq*z'),
                    ]
)

def f(psq,msq):
    # the exact analytical result of the diagram
    msq = complex(msq)
    return 1/psq*(np.log(-psq/msq)+np.log(1-msq/psq))

def fapprox(psq,msq,order):
    # the numerical value of the taylor series to some order of the exact result
    msq = complex(msq)
    s = 0
    for j in range(1,1+order):
        s -= pow(msq/psq,j)/j
    return (s+np.log(-psq/msq))/psq


psq, msq = 1, 0.1
order = 1

print("psq, msq: {}, {}".format(psq,msq))
print("exact result: {:.5f}".format(f(psq,msq)))
print("approximate result to order {}: {:.5f}".format(order,fapprox(psq,msq,order)))
print()

for name, real_parameters, smallness_parameter, li in (
    ("bubble1L_dotted_z_2fp", ['psq','msq','z'], "z", li_z),
    ("bubble1L_dotted_m_2fp", ['psq','msq'], "msq", li_m)
):
    generators_args = psd.loop_integral.loop_regions(
            name = name,
            loop_integral=li,
            smallness_parameter = smallness_parameter,
            expansion_by_regions_order=order)

    psd.code_writer.sum_package(name, [psd.make_package]*len(generators_args), generators_args, li.regulators,
                    requested_orders = [0],
                    real_parameters = real_parameters,
                    complex_parameters = [])