from pySecDec.amplitude_interface import Configuration, Qmc, CubaIntegrator, ContourDeformation

# define the integral directory generated by `sum_pacakge`
config = Configuration('box1L_expansion_by_regions')

# set options for the integrals
config.set_options(
                      Qmc(transform="korobov3", fitfunction="polysingular"),
                      ContourDeformation(),
                  )

# write the config_name.hpp file
config.configure()