from pySecDec.amplitude_interface import Configuration, Qmc, CubaIntegrator, ContourDeformation

config = Configuration("Hjet_g7")

# set options for the integrals
config.set_options(
                    Qmc(transform="korobov3", fitfunction="polysingular"),
                    ContourDeformation(),
                )

# write the config_name.hpp file
config.configure()