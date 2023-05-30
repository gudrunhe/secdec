#include "catch_amalgamated.hpp"
using Catch::Approx;

#include "../secdecutil/coefficient_parser.hpp" // secdecutil::exparse::read_coefficient
#include "../secdecutil/series.hpp" // secdecutil::Series

#include <complex> // std::complex
#include <sstream> // std::istringstream
#include <string> // std::string
#include <vector> // std::vector
#include <fstream> // std::ofstream
#include <cstdio> // std::remove

template<typename T> using nested_series_1_t = secdecutil::Series<T>;
template<typename T> using nested_series_2_t = secdecutil::Series<secdecutil::Series<T>>;

TEST_CASE( "Read a basic coefficient", "[read_coefficient]" ) {

    std::string coeff = "eps^0*(1/2*v3*v4 + 1/2*v3*v4^2*v5*sqrtDelta - 1/2*v3^2*v4^2*sqrtDelta + 1/2*v2*v3^2*v4*sqrtDelta - 1/2*v1*v3*v4*v5*sqrtDelta - 1/2*v1*v2*v3*v4*sqrtDelta)+v1*eps\n";
    std::string coeff_filename = "tmp_basic_coeff.txt";

    std::ofstream coeff_file(coeff_filename);
    coeff_file << coeff;
    coeff_file.close();

    const std::vector<int> required_orders{2};
    const std::vector<std::string> names_of_regulators{"eps"};
    const std::vector<std::string> names_of_real_parameters{"sqrtDelta"};
    const std::vector<std::string> names_of_complex_parameters{"v1","v2","v3","v4","v5"};

    double sqrtDelta = 45;
    std::complex<double> v1 = {-3,1.}, v2 = {-3,1.}, v3 = {-1,0.1}, v4 = {-1,0.1}, v5 = {-1,0.1};
    std::vector<double> real_parameters{sqrtDelta};
    std::vector<std::complex<double>> complex_parameters{v1,v2,v3,v4,v5};

    const nested_series_1_t<std::complex<double>> parsed_coefficient =
        secdecutil::exparse::read_coefficient<nested_series_1_t>
        (
            coeff_filename, required_orders, names_of_regulators,
            names_of_real_parameters, names_of_complex_parameters,
            real_parameters, complex_parameters
        );

    std::remove(coeff_filename.c_str());

    REQUIRE( parsed_coefficient.at(0).real() == Approx(-150.705) );
    REQUIRE( parsed_coefficient.at(0).imag() == Approx(169.55) );
    REQUIRE( parsed_coefficient.at(1).real() == Approx(-3.) );
    REQUIRE( parsed_coefficient.at(1).imag() == Approx(1.) );

    REQUIRE( parsed_coefficient.expansion_parameter == "eps" );

};

TEST_CASE( "Read a multivariate coefficient", "[read_coefficient]" ) {

    std::string coeff = "v1/(v3*eps^3*alp^2 + v2*eps^2*alp^2) + a1/(1 + a2*v2*eps) + (a1 + a2*v1 + a3*eps^2)/(a1*v3 + a2^2*alp + a3^2*eps^3)\n";
    std::string coeff_filename = "tmp_multi_coeff.txt";

    std::ofstream coeff_file(coeff_filename);
    coeff_file << coeff;
    coeff_file.close();

    const std::vector<int> required_orders{3,2};
    const std::vector<std::string> names_of_regulators{"eps","alp"};
    const std::vector<std::string> names_of_real_parameters{"a1","a2","a3"};
    const std::vector<std::string> names_of_complex_parameters{"v1","v2","v3"};

    double a1 = 45.2, a2 = 1.7, a3=-5.7;
    std::complex<double> v1 = {-3.2,1.9}, v2 = {-1.2,3.2}, v3 = {1.1,10.5};
    std::vector<double> real_parameters{a1,a2,a3};
    std::vector<std::complex<double>> complex_parameters{v1,v2,v3};

    const nested_series_2_t<std::complex<double>> parsed_coefficient =
        secdecutil::exparse::read_coefficient<nested_series_2_t>
        (
            coeff_filename, required_orders, names_of_regulators,
            names_of_real_parameters, names_of_complex_parameters,
            real_parameters, complex_parameters
        );

    std::remove(coeff_filename.c_str());

    REQUIRE( parsed_coefficient.at(-2).at(-2).real() == Approx(0.84931506849315068493) );
    REQUIRE( parsed_coefficient.at(-2).at(-2).imag() == Approx(0.68150684931506849315) );
    REQUIRE( parsed_coefficient.at(-1).at(-2).real() == Approx(-3.2878237005066616626) );
    REQUIRE( parsed_coefficient.at(-1).at(-2).imag() == Approx(-0.71130840683054982173) );
    REQUIRE( parsed_coefficient.at(0).at(-2).real() == Approx(10.068256898156121712) );
    REQUIRE( parsed_coefficient.at(0).at(-2).imag() == Approx(-2.5718050239449689859) );
    REQUIRE( parsed_coefficient.at(1).at(-2).real() == Approx(-24.276184562199204521) );
    REQUIRE( parsed_coefficient.at(1).at(-2).imag() == Approx(21.003267754385298019) );
    REQUIRE( parsed_coefficient.at(2).at(-2).real() == Approx(38.10467135848453064) );
    REQUIRE( parsed_coefficient.at(2).at(-2).imag() == Approx(-91.551162521764434669) );
    REQUIRE( parsed_coefficient.at(3).at(-2).real() == Approx(21.043317499911133373) );
    REQUIRE( parsed_coefficient.at(3).at(-2).imag() == Approx(305.60948874155193365) );
    REQUIRE( parsed_coefficient.at(0).at(-1).real() == Approx(0) );
    REQUIRE( parsed_coefficient.at(0).at(-1).imag() == Approx(0) );
    REQUIRE( parsed_coefficient.at(1).at(-1).real() == Approx(0) );
    REQUIRE( parsed_coefficient.at(1).at(-1).imag() == Approx(0) );
    REQUIRE( parsed_coefficient.at(2).at(-1).real() == Approx(0) );
    REQUIRE( parsed_coefficient.at(2).at(-1).imag() == Approx(0) );
    REQUIRE( parsed_coefficient.at(3).at(-1).real() == Approx(0) );
    REQUIRE( parsed_coefficient.at(3).at(-1).imag() == Approx(0) );
    REQUIRE( parsed_coefficient.at(0).at(0).real() == Approx(45.215413085213315146) );
    REQUIRE( parsed_coefficient.at(0).at(0).imag() == Approx(-0.082161107044235084137) );
    REQUIRE( parsed_coefficient.at(1).at(0).real() == Approx(92.208) );
    REQUIRE( parsed_coefficient.at(1).at(0).imag() == Approx(-245.888) );
    REQUIRE( parsed_coefficient.at(2).at(0).real() == Approx(-1149.5276445434609662) );
    REQUIRE( parsed_coefficient.at(2).at(0).imag() == Approx(-1003.2111602669635045) );
    REQUIRE( parsed_coefficient.at(3).at(0).real() == Approx(-7802.5617394423494585) );
    REQUIRE( parsed_coefficient.at(3).at(0).imag() == Approx(4206.8502409315148873) );
    REQUIRE( parsed_coefficient.at(0).at(1).real() == Approx(0.00048514975715804724946) );
    REQUIRE( parsed_coefficient.at(0).at(1).imag() == Approx(0.00014468070415587627843) );
    REQUIRE( parsed_coefficient.at(2).at(1).real() == Approx(-0.000070769087233573255418) );
    REQUIRE( parsed_coefficient.at(2).at(1).imag() == Approx(-0.000014992350652013409759) );
    REQUIRE( parsed_coefficient.at(3).at(1).real() == Approx(-0.000026477162325913921123) );
    REQUIRE( parsed_coefficient.at(3).at(1).imag() == Approx(0.000063650625177592646674) );
    REQUIRE( parsed_coefficient.at(0).at(2).real() == Approx(-1.177577702706851832e-6) );
    REQUIRE( parsed_coefficient.at(0).at(2).imag() == Approx(2.8308757581293128484e-6) );
    REQUIRE( parsed_coefficient.at(2).at(2).real() == Approx(1.3495809617222026119e-7) );
    REQUIRE( parsed_coefficient.at(2).at(2).imag() == Approx(-4.1679845251442038932e-7) );
    REQUIRE( parsed_coefficient.at(3).at(2).real() == Approx(-5.5001297798209142891e-7) );
    REQUIRE( parsed_coefficient.at(3).at(2).imag() == Approx(-2.9946300874021583209e-7) );

    REQUIRE( parsed_coefficient.expansion_parameter == "eps" );
    REQUIRE( parsed_coefficient.at(0).expansion_parameter == "alp" );

};
