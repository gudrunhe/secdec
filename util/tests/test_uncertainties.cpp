#include "catch.hpp"
#include "../secdecutil/uncertainties.hpp"

#include <complex>
#include <sstream>

using dcmplx = std::complex<double>;

TEST_CASE( "Constructor", "[GaussianUncertainty]" ) {

    double one = 1.;
    double two = 2.;

    REQUIRE_NOTHROW(secdecutil::GaussianUncertainty<double>(one,two));

};

TEST_CASE( "Access Fields", "[GaussianUncertainty]" ) {

    double one = 1.;
    double two = 2.;
    auto gu1 = secdecutil::GaussianUncertainty<double>(one,two);

    REQUIRE( gu1.value == one );
    REQUIRE( gu1.uncertainty == two );

};

TEST_CASE( "Unary Operators real", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<double>(3.,5.);

    SECTION( "+" )
    {

        REQUIRE( (+gu1).value == 3. );
        REQUIRE( (+gu1).uncertainty == 5. );

    }

    SECTION( "-" )
    {

        REQUIRE( (-gu1).value == -3. );
        REQUIRE( (-gu1).uncertainty == 5. );

    }

};

TEST_CASE( "Unary Operators complex", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<dcmplx>({3.,1.},{5.,2.});

    SECTION( "+" )
    {

        REQUIRE( (+gu1).value == dcmplx(3.,1.) );
        REQUIRE( (+gu1).uncertainty == dcmplx(5.,2.) );

    }

    SECTION( "-" )
    {

        REQUIRE( (-gu1).value == dcmplx(-3.,-1.) );
        REQUIRE( (-gu1).uncertainty == dcmplx(5.,2.) );

    }

};

TEST_CASE( "Compound Assignment Operators real", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<double>(-1.,2.);
    auto gu2 = secdecutil::GaussianUncertainty<double>(3.,5.);

    SECTION( "+=" )
    {

        gu1 += gu2;
        REQUIRE( gu1.value == Approx(2.) );
        REQUIRE( gu1.uncertainty == Approx( std::sqrt(29.) ) );

    };

    SECTION( "-=" )
    {

        gu1 -= gu2;
        REQUIRE( gu1.value == Approx(-4.) );
        REQUIRE( gu1.uncertainty == Approx( std::sqrt(29.) ));

    };

    SECTION( "*=" )
    {

        gu1 *= gu2;
        REQUIRE( gu1.value == Approx(-3.) );
        REQUIRE( gu1.uncertainty == Approx( 3. * std::sqrt(61./9.) ) );

    };

    SECTION( "/=" )
    {

        gu1 /= gu2;
        REQUIRE( gu1.value == Approx( -1./3.) );
        REQUIRE( gu1.uncertainty == Approx( 1./3. * std::sqrt(61./9.) ) );

    };

};

TEST_CASE( "Compound Assignment Operators complex", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<dcmplx>({-1.,2.},{4.,4.});
    auto gu2 = secdecutil::GaussianUncertainty<dcmplx>({3.,-2.},{3.,3.});

    SECTION( "+=" )
    {

        gu1 += gu2;
        REQUIRE( gu1.value.real() == Approx(2.) );
        REQUIRE( gu1.value.imag() == Approx(0.) );
        REQUIRE( gu1.uncertainty.real() == Approx(5.) );
        REQUIRE( gu1.uncertainty.imag() == Approx(5.) );

    };

    SECTION( "-=" )
    {

        gu1 -= gu2;
        REQUIRE( gu1.value.real() == Approx(-4.) );
        REQUIRE( gu1.value.imag() == Approx( 4.) );
        REQUIRE( gu1.uncertainty.real() == Approx(5.));
        REQUIRE( gu1.uncertainty.imag() == Approx(5.));

    };

    SECTION( "*=" )
    {

        gu1 *= gu2;
        REQUIRE( gu1.value.real() == Approx(1.) );
        REQUIRE( gu1.value.imag() == Approx(8.) );
        REQUIRE( gu1.uncertainty.real() == Approx( std::sqrt( (16. + 1.) * 9.  +  (4. + 9./4.) * 16.) ) );
        REQUIRE( gu1.uncertainty.imag() == Approx( std::sqrt( (16. + 9./4.) * 4.  +  (4. + 1.) * 36.) ) );

    };

    SECTION( "/=" )
    {

        gu1 /= gu2;
        REQUIRE( gu1.value.real() == Approx(-7./13.) );
        REQUIRE( gu1.value.imag() == Approx( 4./13.) );
        REQUIRE( gu1.uncertainty.real() == Approx(   7./13. * std::sqrt( ((16. + 1.) * 9.  +  (4. + 9./4.) * 16.) / 49. + (9.*9.*2. + 9.*4.*2.) / 169.) )   );
        REQUIRE( gu1.uncertainty.imag() == Approx(   4./13. * std::sqrt( ((16. + 9./4.) * 4.  +  (4. + 1.) * 36.) / 16. + (9.*9.*2. + 9.*4.*2.) / 169.) )   );

    };

};

TEST_CASE( "Binary Operators real", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<double>(-6.,8.);
    auto gu2 = secdecutil::GaussianUncertainty<double>(4.,1.);

    SECTION( "+" )
    {

        auto gu3 = gu1 + gu2;
        REQUIRE( gu3.value == Approx(-2.) );
        REQUIRE( gu3.uncertainty == Approx( std::sqrt(65.) ) );

    };

    SECTION( "-" )
    {

        auto gu3 = gu1 - gu2;
        REQUIRE( gu3.value == Approx(-10.) );
        REQUIRE( gu3.uncertainty == Approx( std::sqrt(65.) ));

    };

    SECTION( "*" )
    {

        auto gu3 = gu1 * gu2;
        REQUIRE( gu3.value == Approx(-24.) );
        REQUIRE( gu3.uncertainty == Approx( 24. * std::sqrt(265./144.) ) );

    };

    SECTION( "/" )
    {
        auto gu3 = gu1 / gu2;
        REQUIRE( gu3.value == Approx(-6./4.) );
        REQUIRE( gu3.uncertainty == Approx( 6./4. * std::sqrt(265./144.) ) );

    };
};

TEST_CASE( "Binary Operators complex", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<dcmplx>({-1.,2.},{4.,6.});
    auto gu2 = secdecutil::GaussianUncertainty<dcmplx>({3.,-2.},{3.,8.});

    SECTION( "+" )
    {

        auto gu3 = gu1 + gu2;
        REQUIRE( gu3.value.real() == Approx(2.) );
        REQUIRE( gu3.value.imag() == Approx(0.) );
        REQUIRE( gu3.uncertainty.real() == Approx( 5.) );
        REQUIRE( gu3.uncertainty.imag() == Approx(10.) );

    };

    SECTION( "-" )
    {

        auto gu3 = gu1 - gu2;
        REQUIRE( gu3.value.real() == Approx(-4.) );
        REQUIRE( gu3.value.imag() == Approx( 4.) );
        REQUIRE( gu3.uncertainty.real() == Approx( 5.));
        REQUIRE( gu3.uncertainty.imag() == Approx(10.));

    };

    SECTION( "*" )
    {

        auto gu3 = gu1 * gu2;
        REQUIRE( gu3.value.real() == Approx(1.) );
        REQUIRE( gu3.value.imag() == Approx(8.) );
        REQUIRE( gu3.uncertainty.real() == Approx( std::sqrt( (16. + 1.) * 9.  +  (9. + 16.) * 16.) ) );
        REQUIRE( gu3.uncertainty.imag() == Approx( std::sqrt( (16. + 16.) * 4.  +  (9. + 1.) * 36.) ) );

    };

    SECTION( "/" )
    {
        auto gu3 = gu1 / gu2;
        REQUIRE( gu3.value.real() == Approx(-7./13.) );
        REQUIRE( gu3.value.imag() == Approx( 4./13.) );
        REQUIRE(   gu3.uncertainty.real() == Approx( 7./13. * std::sqrt( ((16. + 1.) * 9.  +  (9. + 16.) * 16.) / 49. + (9.*9.*2. + 4.*4.*4.*4.*2.) / 169.) )   );
        REQUIRE(   gu3.uncertainty.imag() == Approx( 4./13. * std::sqrt( ((16. + 16.) * 4.  +  (9. + 1.) * 36.) / 16. + (9.*9.*2. + 4.*4.*4.*4.*2.) / 169.) )   );

    };

};

TEST_CASE( "Printing real", "[GaussianUncertainty]" ) {

    auto gu = secdecutil::GaussianUncertainty<double>(1.0,0.5);
    std::stringstream output;
    output << gu;
    REQUIRE( output.str() == "1 +/- 0.5" );

};

TEST_CASE( "Printing complex", "[GaussianUncertainty]" ) {

    auto gu = secdecutil::GaussianUncertainty<std::complex<double>>({1.0,1.5},{0.5,0.1});
    std::stringstream output;
    output << gu;
    REQUIRE( output.str() == "(1,1.5) +/- (0.5,0.1)" );

};
