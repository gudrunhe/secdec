#include "%(name)s.hpp"
#include <vector> 
#include <secdecutil/uncertainties.hpp> // secdecutil::UncorrelatedDeviation
#include <string>
#include <sstream>
#define INTEGRAL_NAME %(name)s

/*
 * First perform replacements, then put item into quotes.
 * With only a single macro, the preprocessor does not
 * consider replacements in "item".
 */
#define EXPAND_STRINGIFY(item) STRINGIFY(item)
#define STRINGIFY(item) #item

using namespace INTEGRAL_NAME;

template<typename T> using amplitudes_t = std::vector<%(name)s::nested_series_t<T>>;

extern "C"
{
    std::string * allocate_string()
    {
        return new std::string;
    }

    const char * string2charptr(std::string * str)
    {
        return str->c_str();
    }

    void free_string(std::string * strptr)
    {
        delete strptr;
    }

    void compute_integral
    (
        std::string * integral_result_strptr,
        const double real_parameters_input[], // real parameters
        const double complex_parameters_input[], // complex parameters serialized as real(x0), imag(x0), real(x1), imag(x1), ...
        int verbose = 1
    )
    {
        size_t i;
        
        // read real parameters
        std::vector<%(name)s::real_t> real_parameters(number_of_real_parameters);
        for (i=0 ; i<number_of_real_parameters ; ++i)
            real_parameters[i] = real_parameters_input[i];

        // read complex parameters
        std::vector<%(name)s::complex_t> complex_parameters(number_of_complex_parameters);
        for (i=0 ; i<number_of_complex_parameters ; ++i)
            complex_parameters[i] = complex_t(complex_parameters_input[2*i],complex_parameters_input[2*i + 1]);
        
        // Construct the amplitudes. Further options for the individual integrals can be set in the
        // corresponding "<integral_name>_weighted_integral.cpp" file in the "src/" directory
        std::vector<%(name)s::nested_series_t<%(name)s::sum_t>> unwrapped_amplitudes =
            %(name)s::make_amplitudes(real_parameters, complex_parameters);

        // pack amplitude into handler
        %(name)s::handler_t<amplitudes_t> amplitudes
        (
            unwrapped_amplitudes
            // further optional arguments: epsrel, epsabs, maxeval, mineval, maxincreasefac, min_epsrel, min_epsabs, max_epsrel, max_epsabs
        );
        amplitudes.verbose = verbose;

        // compute the amplitude
        const std::vector<%(name)s::nested_series_t<secdecutil::UncorrelatedDeviation<%(name)s::integrand_return_t>>> result = amplitudes.evaluate();

        std::stringstream sstream;

        sstream.precision(std::numeric_limits<%(name)s::real_t>::max_digits10); // force enough digits to ensure unique recreation
        sstream << std::scientific; // stringify floats as #.#e#

        sstream.str("");
        for(const auto amp : result)
            sstream << amp;

        *integral_result_strptr = sstream.str();
    }

    void get_integral_info(std::string * str_ptr)
    {
        std::stringstream sstream;

        // fix output formatting
        sstream.precision(std::numeric_limits<real_t>::max_digits10); // force enough digits to ensure unique recreation
        sstream << std::scientific; // stringify floats as #.#e#

        sstream << "name = " << EXPAND_STRINGIFY(INTEGRAL_NAME) << std::endl;

        sstream << "number_of_regulators = " << number_of_regulators << std::endl;
        sstream << "names_of_regulators =";
        for ( const auto& name : names_of_regulators )
            sstream << " " << name;
        sstream << std::endl;

        sstream << "number_of_real_parameters = " << number_of_real_parameters << std::endl;
        sstream << "names_of_real_parameters =";
        for ( const auto& name : names_of_real_parameters )
            sstream << " " << name;
        sstream << std::endl;

        sstream << "number_of_complex_parameters = " << number_of_complex_parameters << std::endl;
        sstream << "names_of_complex_parameters =";
        for ( const auto& name : names_of_complex_parameters )
            sstream << " " << name;
        // sstream << std::endl;

        *str_ptr = sstream.str();

    }
}
