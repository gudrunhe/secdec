#ifndef %(name)s_config_%(name)s_hpp_included
#define %(name)s_config_%(name)s_hpp_included

/*
 * Set contour-deformation and integrator properties. These can be set per
 * integral by defining an Options struct in the appropriate namespace
 * (see below). Integranls for which no special Options struct is defined
 * use the Options struct in the parent namespace of all the integrands.
 */

namespace %(name)s
{
    // This struct contains the options for all integrals which do not set idividual options.
    struct Options
    {
        // whether the integrator to be used supports CUDA (out of the builtin integrators, only the Qmc does)
        constexpr static bool cuda_compliant_integrator = false; // TODO: generate with python script

        // integral settings
        template<typename integrand_t, unsigned long long int maximal_number_of_integration_variables>
        struct Integral
        {
            // the integrator (Cuhre, Divonne, Qmc, Suave, Vegas)
            using integrator_t = secdecutil::cuba::Vegas<integrand_return_t>; // TODO: generate with python script

            // instantiation and options of the integrator
            static std::shared_ptr<integrator_t> configure_integrator()
            {
                std::shared_ptr<integrator_t> integrator = std::make_shared<integrator_t>();
                // further options can be set here, e.g. generating vectors for Qmc etc
                //integrator->nincrease = 1e5;
                // TODO: generate with python script
                return integrator;
            }

            // the kind of integral (CubaIntegral, QmcIntegral, <CustomIntegral>)
            using integral_t = secdecutil::amplitude::CubaIntegral<integrand_return_t,real_t,integrator_t,integrand_t>; // TODO: generate with python script
        };

        // parameters for the contour deformation
        struct ContourDeformation
        {
            constexpr static unsigned number_of_presamples = 10000; // TODO: generate with python script
            constexpr static real_t deformation_parameters_maximum = 1.; // TODO: generate with python script
            constexpr static real_t deformation_parameters_minimum = 1.e-5; // TODO: generate with python script
            constexpr static real_t deformation_parameters_decrease_factor = 0.9; // TODO: generate with python script
        };
    };


    // TODO: generate integral-specific options with python script
    /*
    // The global options above can be overruled for a particular integral
    // by defining an Options struct in the appropriate namespace.
    namespace SUB_INTEGRAL_NAME
    {
        struct Options
        {
            ...
        }
    };
    */
}

#endif
