#include "catch.hpp"
#include "../secdecutil/uncertainties.hpp"

#include <cmath>
#include <iostream>

TEST_CASE( "Constructor", "[GaussianUncertainty]" ) {

    double one = 1.;
    double two = 2.;

    REQUIRE_NOTHROW(secdecutil::GaussianUncertainty<double>(one,two));

};

TEST_CASE( "Access Fields", "[GaussianUncertainty]" ) {

    double one = 1.;
    double two = 2.;
    auto gu1 = secdecutil::GaussianUncertainty<double>(one,two);

    REQUIRE( gu1.value == Approx(one) );
    REQUIRE( gu1.uncertainty == Approx(two) );

};

TEST_CASE( "Unary Operators", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<double>(3.,5.);

    SECTION( "+" )
    {

        REQUIRE( gu1.value == Approx(3.) );
        REQUIRE( gu1.uncertainty == Approx(5.) );

    }

    SECTION( "-" )
    {

        REQUIRE( (-gu1).value == Approx(-3.) );
        REQUIRE( (-gu1).uncertainty == Approx(5.) );

    }

};

TEST_CASE( "Compound Assignment Operators", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<double>(-1.,2.);
    auto gu2 = secdecutil::GaussianUncertainty<double>(3.,5.);

    SECTION( "+=" )
    {

        gu1 += gu2;
        REQUIRE( gu1.value == Approx(2.) );
        REQUIRE( gu1.uncertainty == Approx( sqrt(29.) ) );

    };

    SECTION( "-=" )
    {

        gu1 -= gu2;
        REQUIRE( gu1.value == Approx(-4.) );
        REQUIRE( gu1.uncertainty == Approx( sqrt(29.) ));

    };

    SECTION( "*=" )
    {

        gu1 *= gu2;
        REQUIRE( gu1.value == Approx(-3.) );
        REQUIRE( gu1.uncertainty == Approx( 3. * sqrt(61./9.) ) );
        
    };

    SECTION( "/=" )
    {

        gu1 /= gu2;
        REQUIRE( gu1.value == Approx( -1./3.) );
        REQUIRE( gu1.uncertainty == Approx( 1./3. * sqrt(61./9.) ) );

    };

};

TEST_CASE( "Binary Operators", "[GaussianUncertainty]" ) {

    auto gu1 = secdecutil::GaussianUncertainty<double>(-6.,8.);
    auto gu2 = secdecutil::GaussianUncertainty<double>(4.,1.);

    SECTION( "+" )
    {

        auto gu3 = gu1 + gu2;
        REQUIRE( gu3.value == Approx(-2.) );
        REQUIRE( gu3.uncertainty == Approx( sqrt(65.) ) );

    };

    SECTION( "-" )
    {

        auto gu3 = gu1 - gu2;
        REQUIRE( gu3.value == Approx(-10.) );
        REQUIRE( gu3.uncertainty == Approx( sqrt(65.) ));

    };

    SECTION( "*" )
    {

        auto gu3 = gu1 * gu2;
        REQUIRE( gu3.value == Approx(-24.) );
        REQUIRE( gu3.uncertainty == Approx( 24. * sqrt(265./144.) ) );

    };

    SECTION( "/" )
    {
        auto gu3 = gu1 / gu2;
        REQUIRE( gu3.value == Approx(-6./4.) );
        REQUIRE( gu3.uncertainty == Approx( 6./4. * sqrt(265./144.) ) );
        
    };
    
};
