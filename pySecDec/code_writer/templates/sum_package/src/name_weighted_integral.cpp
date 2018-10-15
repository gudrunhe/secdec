#include <cassert> // assert
#include <fstream> // std::ifstream
#include <memory> // std::shared_ptr, std::make_shared
#include <numeric> // std::accumulate
#include <sstream> // std::stringstream
#include <string> // std::string
#include <vector> // std::vector

#include <secdecutil/amplitude.hpp> // secdecutil::amplitude::Integral, secdecutil::amplitude::CubaIntegral, secdecutil::amplitude::QmcIntegral
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply
#include <secdecutil/ginac_coefficient_parser.hpp> // secdecutil::ginac::read_coefficient
#include <secdecutil/integrators/cuba.hpp> // secdecutil::cuba::Vegas, secdecutil::cuba::Suave, secdecutil::cuba::Cuhre, secdecutil::cuba::Divonne
#include <secdecutil/integrators/qmc.hpp> // secdecutil::integrators::Qmc
#include <secdecutil/series.hpp> // secdecutil::Series

#include <ginac/ginac.h>

#define QUOTE(ARG) #ARG
#define QUOTE_EXPAND(ARG) QUOTE(ARG)
#define CAT(ARG1,ARG2) ARG1##ARG2
#define CAT_EXPAND(ARG1,ARG2) CAT(ARG1,ARG2)

#define MAIN_INTEGRAL_NAME %(name)s
#define SUB_INTEGRAL_NAME %(sub_integral_name)s

#include QUOTE_EXPAND(MAIN_INTEGRAL_NAME.hpp)
#include QUOTE_EXPAND(SUB_INTEGRAL_NAME/SUB_INTEGRAL_NAME.hpp)
#include QUOTE_EXPAND(CAT_EXPAND(SUB_INTEGRAL_NAME,_weighted_integral.hpp))

#define CONTOURDEF_PARAMETERS \
    // put parameters for contour deformation here: number_of_presamples, deformation_parameters_maximum, deformation_parameters_minimum, deformation_parameters_decrease_factor

#define CONFIGURE_INTEGRATOR \
    using integrator_t = secdecutil::integrators::Qmc< \
                                                         integrand_return_t, \
                                                         maximal_number_of_integration_variables, \
                                                         ::integrators::transforms::Korobov<3>::type, \
                                                         integrand_t, \
                                                         ::integrators::fitfunctions::PolySingular::type \
                                                     >; \
    auto integrator = std::make_shared<integrator_t>(); \
    // further options can be set here, e.g.: integrator->verbose = true;

#define MAKE_WEIGHTED_INTEGRAL_BODY(USING_INTEGRAND_T,CUDA) \
    using namespace SUB_INTEGRAL_NAME; \
    USING_INTEGRAND_T; \
    using integral_t = secdecutil::amplitude::Integral<integrand_return_t,real_t>; \
 \
    CONFIGURE_INTEGRATOR \
 \
    std::vector<nested_series_t<integrand_t>> raw_integrands = \
        make_##CUDA##integrands \
        ( \
            real_parameters, \
            complex_parameters \
            CONTOURDEF_PARAMETERS \
        ); \
 \
    const std::function<sum_t(const integrand_t& sector)> convert_integrands = \
        [ integrator ] (const integrand_t& sector) -> sum_t \
        { \
            return { /* constructor of std::vector */ \
                        { /* constructor of WeightedIntegral */  \
                            std::make_shared< \
                                                secdecutil::amplitude::QmcIntegral< \
                                                                                       integrand_return_t, \
                                                                                       real_t, \
                                                                                       integrator_t, \
                                                                                       integrand_t \
                                                                                  > \
                                            >(integrator,sector) \
                        } \
                   }; \
        }; \
 \
    const std::vector<nested_series_t<sum_t>> integrals = deep_apply(raw_integrands, convert_integrands); \
    nested_series_t<sum_t> amplitude = std::accumulate(++integrals.begin(), integrals.end(), *integrals.begin() ); \
    amplitude *= prefactor(real_parameters,complex_parameters) * coefficient(real_parameters,complex_parameters); \
    return amplitude;

namespace SUB_INTEGRAL_NAME
{
    const std::vector<int> lowest_coefficient_orders{%(lowest_coefficient_orders)s};

    static std::vector<int> compute_required_orders()
    {
        assert(MAIN_INTEGRAL_NAME::requested_orders.size() == SUB_INTEGRAL_NAME::requested_orders.size());
        size_t number_of_regulators = MAIN_INTEGRAL_NAME::requested_orders.size();
        std::vector<int> orders; orders.reserve( MAIN_INTEGRAL_NAME::requested_orders.size() );
        for(size_t i = 0; i < number_of_regulators; ++i)
            orders.push_back(
                                  MAIN_INTEGRAL_NAME::requested_orders.at(i) + 1
                                - SUB_INTEGRAL_NAME::lowest_coefficient_orders.at(i)
                                - SUB_INTEGRAL_NAME::lowest_orders.at(i)
                                - SUB_INTEGRAL_NAME::lowest_prefactor_orders.at(i)
                            );
        return orders;
    }

    nested_series_t<complex_t> coefficient(const std::vector<real_t>& real_parameters, const std::vector<complex_t>& complex_parameters)
    {
        std::ifstream coeffile(QUOTE_EXPAND(lib/CAT_EXPAND(SUB_INTEGRAL_NAME,_coefficient.txt)));
        assert( coeffile.is_open() );
        return secdecutil::ginac::read_coefficient<nested_series_t>
               (
                    coeffile, compute_required_orders(),
                    names_of_regulators, names_of_real_parameters, names_of_complex_parameters,
                    real_parameters, complex_parameters
               );
    }
};

namespace MAIN_INTEGRAL_NAME
{
    // convert return value of "make_integrands", "prefactor", and "coefficient" to "nested_series_t<sum_t>"
    nested_series_t<sum_t> CAT_EXPAND(make_weighted_integral_,SUB_INTEGRAL_NAME)
    (
        const std::vector<real_t>& real_parameters,
        const std::vector<complex_t>& complex_parameters
    )
    {
        #ifdef SECDEC_WITH_CUDA
            MAKE_WEIGHTED_INTEGRAL_BODY(using integrand_t = SUB_INTEGRAL_NAME::cuda_integrand_t, cuda_)
        #else
            MAKE_WEIGHTED_INTEGRAL_BODY(using SUB_INTEGRAL_NAME::integrand_t,)
        #endif
    }
};
