#include "catch.hpp"
#include "../secdecutil/uncertainties.hpp"

#include <complex>
#include <sstream>

#ifdef SECDEC_WITH_CUDA
    #include <thrust/complex.h>
    using dcmplx = thrust::complex<double>;
#else
    using dcmplx = std::complex<double>;
#endif

TEST_CASE( "Access Fields", "[UncorrelatedDeviation]" ) {

    double one = 1.;
    double two = 2.;

    SECTION( "with uncertainty" )
    {

        auto gu1 = secdecutil::UncorrelatedDeviation<double>(one,two);

        REQUIRE( gu1.value == one );
        REQUIRE( gu1.uncertainty == two );

    }

    SECTION( "zero uncertainty" )
    {

        auto gu1 = secdecutil::UncorrelatedDeviation<double>(one);

        REQUIRE( gu1.value == one );
        REQUIRE( gu1.uncertainty == 0 );

    }


};

TEST_CASE( "Unary Operators real", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<double>(3.,5.);

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

TEST_CASE( "Unary Operators complex", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<dcmplx>({3.,1.},{5.,2.});

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

TEST_CASE( "Compound Assignment Operators real", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<double>(-1.,2.);
    auto gu2 = secdecutil::UncorrelatedDeviation<double>(3.,5.);

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
        REQUIRE( gu1.uncertainty == Approx( std::sqrt(161.) ) );

    };

};

TEST_CASE( "Compound Assignment Operators real zero uncertainty", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<double>(-1.,2.);
    auto gu2 = secdecutil::UncorrelatedDeviation<double>(3.);

    SECTION( "*=" )
    {

        gu1 *= gu2;
        REQUIRE( gu1.value == Approx(-3.) );
        REQUIRE( gu1.uncertainty == Approx( 6. ) );

    };

};

TEST_CASE( "Compound Assignment Operators real plain number", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<double>(-1.,2.);

    SECTION( "+=" )
    {

        gu1 += 3.;
        REQUIRE( gu1.value == Approx(2.) );
        REQUIRE( gu1.uncertainty == Approx( 2. ) );

    };

    SECTION( "+=" )
    {

        gu1 -= 3.;
        REQUIRE( gu1.value == Approx(-4.) );
        REQUIRE( gu1.uncertainty == Approx( 2. ) );

    };

    SECTION( "*=" )
    {

        gu1 *= 3.;
        REQUIRE( gu1.value == Approx(-3.) );
        REQUIRE( gu1.uncertainty == Approx( 6. ) );

    };

    SECTION( "/=" )
    {

        gu1 /= 3.;
        REQUIRE( gu1.value == Approx( -1./3.) );
        REQUIRE( gu1.uncertainty == Approx( 2./3. ) );

    };

};

TEST_CASE( "Compound Assignment Operators complex", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<dcmplx>({-1.,2.},{4.,4.});
    auto gu2 = secdecutil::UncorrelatedDeviation<dcmplx>({3.,-2.},{3.,3.});

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
        REQUIRE( gu1.uncertainty.real() == Approx(5.) );
        REQUIRE( gu1.uncertainty.imag() == Approx(5.) );

    };

    SECTION( "*=" )
    {

        gu1 *= gu2;
        REQUIRE( gu1.value.real() == Approx(1.) );
        REQUIRE( gu1.value.imag() == Approx(8.) );
        REQUIRE( gu1.uncertainty.real() == Approx( std::sqrt(541.) ) );
        REQUIRE( gu1.uncertainty.imag() == Approx( std::sqrt(541.) ) );

    };

};

TEST_CASE( "Compound Assignment Operators complex plain number", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<dcmplx>({-1.,2.},{4.,4.});

    SECTION( "+=" )
    {

        gu1 += dcmplx{3.,-2.};
        REQUIRE( gu1.value.real() == Approx(2.) );
        REQUIRE( gu1.value.imag() == Approx(0.) );
        REQUIRE( gu1.uncertainty.real() == Approx(4.) );
        REQUIRE( gu1.uncertainty.imag() == Approx(4.) );

    };

    SECTION( "-=" )
    {

        gu1 -= dcmplx{3.,-2.};
        REQUIRE( gu1.value.real() == Approx(-4.) );
        REQUIRE( gu1.value.imag() == Approx( 4.) );
        REQUIRE( gu1.uncertainty.real() == Approx(4.));
        REQUIRE( gu1.uncertainty.imag() == Approx(4.));

    };

    SECTION( "*=" )
    {

        gu1 *= dcmplx{3.,-2.};
        REQUIRE( gu1.value.real() == Approx(1.) );
        REQUIRE( gu1.value.imag() == Approx(8.) );
        REQUIRE( gu1.uncertainty.real() == Approx( std::sqrt( (16. + 0.) * 9.  +  (4. + 0./4.) * 16.) ) );
        REQUIRE( gu1.uncertainty.imag() == Approx( std::sqrt( (16. + 0./4.) * 4.  +  (4. + 0.) * 36.) ) );

    };

    SECTION( "/=" )
    {

        gu1 /= dcmplx{3.,-2.};
        REQUIRE( gu1.value.real() == Approx(-7./13.) );
        REQUIRE( gu1.value.imag() == Approx( 4./13.) );
        REQUIRE( gu1.uncertainty.real() == Approx(   7./13. * std::sqrt( ((16. + 0.) * 9.  +  (4. + 0./4.) * 16.) / 49. + (0.*9.*2. + 0.*4.*2.) / 169.) )   );
        REQUIRE( gu1.uncertainty.imag() == Approx(   4./13. * std::sqrt( ((16. + 0./4.) * 4.  +  (4. + 0.) * 36.) / 16. + (0.*9.*2. + 0.*4.*2.) / 169.) )   );

    };

};

TEST_CASE( "Binary Operators real", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<double>(-6.,8.);
    auto gu2 = secdecutil::UncorrelatedDeviation<double>(4.,1.);

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
        REQUIRE( gu3.uncertainty == Approx( 2. * std::sqrt(281.) ) );

    };

};

TEST_CASE( "Binary Operators real plain number", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<double>(-6.,8.);

    SECTION( "+" )
    {

        auto gu2 = 4. + gu1;
        REQUIRE( gu2.value == Approx(-2.) );
        REQUIRE( gu2.uncertainty == Approx( 8. ) );

        auto gu3 = gu1 + 4.;
        REQUIRE( gu3.value == Approx(-2.) );
        REQUIRE( gu3.uncertainty == Approx( 8. ) );

    };

    SECTION( "-" )
    {

        auto gu2 = 4. - gu1;
        REQUIRE( gu2.value == Approx(+10.) );
        REQUIRE( gu2.uncertainty == Approx( 8. ));

        auto gu3 = gu1 - 4.;
        REQUIRE( gu3.value == Approx(-10.) );
        REQUIRE( gu3.uncertainty == Approx( 8. ));

    };

    SECTION( "*" )
    {

        auto gu2 = 4. * gu1;
        REQUIRE( gu2.value == Approx(-24.) );
        REQUIRE( gu2.uncertainty == Approx( 24. * std::sqrt( (8.*8.) / (6.*6.) ) ) );

        auto gu3 = gu1 * 4.;
        REQUIRE( gu3.value == Approx(-24.) );
        REQUIRE( gu3.uncertainty == Approx( 24. * std::sqrt( (8.*8.) / (6.*6.) ) ) );

    };

    SECTION( "/" )
    {

        // division by UncorrelatedDeviation is not defined

        auto gu3 = gu1 / 4.;
        REQUIRE( gu3.value == Approx(-6./4.) );
        REQUIRE( gu3.uncertainty == Approx( 6./4. * std::sqrt( (8.*8.) / (6.*6.) ) ) );

    };
};

TEST_CASE( "Implicit conversion", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<std::complex<double>>(-6.,8.);
    auto gu2 = secdecutil::UncorrelatedDeviation<double>(4);

    SECTION( "+" )
    {

        auto gu3 = gu1 + gu2;
        REQUIRE( gu3.value.real() == Approx(-2.) );
        REQUIRE( gu3.value.imag() == Approx( 0.) );
        REQUIRE( gu3.uncertainty.real() == Approx( 8. ) );
        REQUIRE( gu3.uncertainty.imag() == Approx( 0. ) );

        auto gu4 = gu2 + gu1;
        REQUIRE( gu4.value.real() == Approx(-2.) );
        REQUIRE( gu4.value.imag() == Approx( 0.) );
        REQUIRE( gu4.uncertainty.real() == Approx( 8. ) );
        REQUIRE( gu4.uncertainty.imag() == Approx( 0. ) );

    };

    SECTION( "-" )
    {

        auto gu3 = gu2 - gu1;
        REQUIRE( gu3.value.real() == Approx(10.) );
        REQUIRE( gu3.value.imag() == Approx( 0.) );
        REQUIRE( gu3.uncertainty.real() == Approx( 8. ) );
        REQUIRE( gu3.uncertainty.imag() == Approx( 0. ) );

        auto gu4 = gu1 - gu2;
        REQUIRE( gu4.value.real() == Approx(-10.) );
        REQUIRE( gu4.value.imag() == Approx(  0.) );
        REQUIRE( gu4.uncertainty.real() == Approx( 8. ) );
        REQUIRE( gu4.uncertainty.imag() == Approx( 0. ) );

    };

    SECTION( "*" )
    {

        auto gu3 = gu1 * gu2;
        REQUIRE( gu3.value.real() == Approx(-24.) );
        REQUIRE( gu3.value.imag() == Approx(  0.) );
        REQUIRE( gu3.uncertainty.real() == Approx( 24. * std::sqrt( (8.*8.) / (6.*6.) ) ) );
        REQUIRE( gu3.uncertainty.imag() == Approx( 0. ) );

        auto gu4 = gu2 * gu1;
        REQUIRE( gu4.value.real() == Approx(-24.) );
        REQUIRE( gu4.value.imag() == Approx(  0.) );
        REQUIRE( gu4.uncertainty.real() == Approx( 24. * std::sqrt( (8.*8.) / (6.*6.) ) ) );
        REQUIRE( gu4.uncertainty.imag() == Approx( 0. ) );

    };

};

TEST_CASE( "Explicit conversion", "[UncorrelatedDeviation]" ) {

    secdecutil::UncorrelatedDeviation<int> real_deviation = 1;

    secdecutil::UncorrelatedDeviation<std::complex<int>> complex_deviation_0 = 1;
    secdecutil::UncorrelatedDeviation<std::complex<int>> complex_deviation_1 = std::complex<int>{1,2};
    auto complex_deviation_2 = secdecutil::UncorrelatedDeviation<std::complex<int>>({1,1},2);
    auto complex_deviation_3 = secdecutil::UncorrelatedDeviation<std::complex<int>>(1,{2,1});

    secdecutil::UncorrelatedDeviation<std::complex<int>> complex_deviation_4 = real_deviation;

    REQUIRE( real_deviation.value                   == 1 );
    REQUIRE( real_deviation.uncertainty             == 0 );

    REQUIRE( complex_deviation_0.value.real()       == 1 );
    REQUIRE( complex_deviation_0.value.imag()       == 0 );
    REQUIRE( complex_deviation_0.uncertainty.real() == 0 );
    REQUIRE( complex_deviation_0.uncertainty.imag() == 0 );

    REQUIRE( complex_deviation_1.value.real()       == 1 );
    REQUIRE( complex_deviation_1.value.imag()       == 2 );
    REQUIRE( complex_deviation_1.uncertainty.real() == 0 );
    REQUIRE( complex_deviation_1.uncertainty.imag() == 0 );

    REQUIRE( complex_deviation_2.value.real()       == 1 );
    REQUIRE( complex_deviation_2.value.imag()       == 1 );
    REQUIRE( complex_deviation_2.uncertainty.real() == 2 );
    REQUIRE( complex_deviation_2.uncertainty.imag() == 0 );

    REQUIRE( complex_deviation_3.value.real()       == 1 );
    REQUIRE( complex_deviation_3.value.imag()       == 0 );
    REQUIRE( complex_deviation_3.uncertainty.real() == 2 );
    REQUIRE( complex_deviation_3.uncertainty.imag() == 1 );

    REQUIRE( complex_deviation_4.value.real()       == 1 );
    REQUIRE( complex_deviation_4.value.imag()       == 0 );
    REQUIRE( complex_deviation_4.uncertainty.real() == 0 );
    REQUIRE( complex_deviation_4.uncertainty.imag() == 0 );

};

TEST_CASE( "Binary Operators complex", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<dcmplx>({-1.,2.},{4.,6.});
    auto gu2 = secdecutil::UncorrelatedDeviation<dcmplx>({3.,-2.},{3.,8.});

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
        REQUIRE( gu3.uncertainty.real() == Approx( std::sqrt(3001.) ) );
        REQUIRE( gu3.uncertainty.imag() == Approx( 6. * std::sqrt(51.) ) );

    };

};

TEST_CASE( "Binary Operators complex plain number", "[UncorrelatedDeviation]" ) {

    auto gu1 = secdecutil::UncorrelatedDeviation<dcmplx>({-1.,2.},{4.,6.});

    SECTION( "+" )
    {

        auto gu3 = gu1 + dcmplx{3.,-2.};
        REQUIRE( gu3.value.real() == Approx(2.) );
        REQUIRE( gu3.value.imag() == Approx(0.) );
        REQUIRE( gu3.uncertainty.real() == Approx(4.) );
        REQUIRE( gu3.uncertainty.imag() == Approx(6.) );

    };

    SECTION( "-" )
    {

        auto gu3 = gu1 - dcmplx{3.,-2.};
        REQUIRE( gu3.value.real() == Approx(-4.) );
        REQUIRE( gu3.value.imag() == Approx( 4.) );
        REQUIRE( gu3.uncertainty.real() == Approx(4.));
        REQUIRE( gu3.uncertainty.imag() == Approx(6.));

    };

    SECTION( "*" )
    {

        auto gu3 = gu1 * dcmplx{3.,-2.};
        REQUIRE( gu3.value.real() == Approx(1.) );
        REQUIRE( gu3.value.imag() == Approx(8.) );
        REQUIRE( gu3.uncertainty.real() == Approx( std::sqrt( (16. + 0.) * 9.  +  (0. + 9.) * 16.) ) );
        REQUIRE( gu3.uncertainty.imag() == Approx( std::sqrt( (16. + 0.) * 4.  +  (0. + 9.) * 36.) ) );

    };

    SECTION( "/" )
    {
        auto gu3 = gu1 / dcmplx{3.,-2.};
        REQUIRE( gu3.value.real() == Approx(-7./13.) );
        REQUIRE( gu3.value.imag() == Approx( 4./13.) );
        REQUIRE(   gu3.uncertainty.real() == Approx( 7./13. * std::sqrt( ((16. + 0.) * 9.  +  (0. + 9.) * 16.) / 49.) )   );
        REQUIRE(   gu3.uncertainty.imag() == Approx( 4./13. * std::sqrt( ((0. + 16.) * 4.  +  (9. + 0.) * 36.) / 16.) )   );

    };

};

TEST_CASE( "Printing real", "[UncorrelatedDeviation]" ) {

    auto gu = secdecutil::UncorrelatedDeviation<double>(1.0,0.5);
    std::stringstream output;
    output << gu;
    REQUIRE( output.str() == "1 +/- 0.5" );

};

TEST_CASE( "Printing complex", "[UncorrelatedDeviation]" ) {

    auto gu = secdecutil::UncorrelatedDeviation<std::complex<double>>({1.0,1.5},{0.5,0.1});
    std::stringstream output;
    output << gu;
    REQUIRE( output.str() == "(1,1.5) +/- (0.5,0.1)" );

};
