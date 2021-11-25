#include "catch.hpp"
#include "../secdecutil/coefficient_parser.hpp" // secdecutil::exparse::read_coefficient
#include "../secdecutil/series.hpp" // secdecutil::Series

#include <complex> // std::complex
#include <sstream> // std::istringstream
#include <string> // std::string
#include <vector> // std::vector

TEST_CASE( "Read a basic coefficient", "[read_coefficient]" ) {

    std::istringstream coeff_stream
    (
        "\n\
        regulator_factor = 1\n\
        *eps^0\n\
        ;\n\
        numerator = 1\n\
        ,1/2*v3*v4 + 1/2*v3*v4^2*v5*sqrtDelta - 1/2*v3^2*v4^2*sqrtDelta + 1/2*\n\
          v2*v3^2*v4*sqrtDelta - 1/2*v1*v3*v4*v5*sqrtDelta - 1/2*v1*v2*v3*v4\n\
          *sqrtDelta\n\
        ;\n\
        denominator = 1\n\
        ,1\n\
        ;\n\
        "
    );


    const std::vector<int> required_orders{2};
    const std::vector<std::string> names_of_regulators{"eps"};
    const std::vector<std::string> names_of_real_parameters{"sqrtDelta"};
    const std::vector<std::string>& names_of_complex_parameters{"v1","v2","v3","v4","v5"};

    double sqrtDelta = 45;
    std::vector<double> real_parameters{sqrtDelta};

    std::complex<double> v1 = {-3,1.}, v2 = {-3,1.}, v3 = {-1,0.1}, v4 = {-1,0.1}, v5 = {-1,0.1};
    std::vector<std::complex<double>> complex_parameters{v1,v2,v3,v4,v5};

    const secdecutil::Series<std::complex<double>> parsed_coefficient =
        secdecutil::exparse::read_coefficient<secdecutil::Series>
        (
            coeff_stream, required_orders, names_of_regulators,
            names_of_real_parameters, names_of_complex_parameters,
            real_parameters, complex_parameters
        );

    REQUIRE( parsed_coefficient.get_order_min() == 0 );
    REQUIRE( parsed_coefficient.get_order_max() == 1 );

    REQUIRE( parsed_coefficient.at(0).real() == Approx(-150.705) );
    REQUIRE( parsed_coefficient.at(0).imag() == Approx(169.55) );
    REQUIRE( parsed_coefficient.at(1).real() == Approx(0) );
    REQUIRE( parsed_coefficient.at(1).imag() == Approx(0) );

    REQUIRE( parsed_coefficient.expansion_parameter == "eps" );

};
