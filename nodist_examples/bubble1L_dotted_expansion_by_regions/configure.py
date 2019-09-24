from pySecDec.amplitude_interface import Configuration, Qmc, CubaIntegrator, ContourDeformation
from os.path import isdir

for name in "bubble1L_dotted_z_3fp", "bubble1L_dotted_m_3fp", \
            "bubble1L_dotted_z_2fp", "bubble1L_dotted_m_2fp":
    if not isdir(name):
        continue
    # give the integral directory generated by `sum_pacakge`
    config = Configuration(name)

    # set options for the integrals
    config.set_options(
                        Qmc(transform="korobov3", fitfunction="polysingular"),
                        ContourDeformation(),
                    )

    # write the config_box2L_pp.hpp file
    config.configure()
