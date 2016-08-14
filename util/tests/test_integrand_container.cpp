#include <numeric>
#include <array>

#include "catch.hpp"
#include "../secdecutil/integrand_container.hpp"

TEST_CASE( "Check Access", "[IntegrandContainer]" ) {

    const std::function<int(int)> func = [] (int i) { return i+2; };
    auto ic = secdecutil::IntegrandContainer<int, int>(3,func);

    SECTION( "Accessing fields" ) {

        REQUIRE( ic.number_of_integration_variables == 3 );
        REQUIRE( ic.integrand(10) == 10+2 );

    };

};

TEST_CASE( "Binary operators + and -", "[IntegrandContainer]" ) {

    std::function<int(int)> func1 = [] (int i) { return i+2; };
    std::function<int(int)> func2 = [] (int i) { return i+5; };

    auto ic1 = secdecutil::IntegrandContainer<int, int>(3,func1);
    auto ic2 = secdecutil::IntegrandContainer<int, int>(4,func2);
    auto ic3 = ic1 + ic2;
    auto ic4 = ic1 - ic2;

    SECTION ( " + " ) {

        REQUIRE( ic3.number_of_integration_variables == 4);
        REQUIRE( ic3.integrand(10) == 10+2+10+5);

    };

    SECTION ( " += " ) {

        ic1 += ic2;

        REQUIRE( ic1.number_of_integration_variables == 4);
        REQUIRE( ic1.integrand(10) == 10+2+10+5);

    };

    SECTION ( " - " ) {

        REQUIRE( ic4.number_of_integration_variables == 4);
        REQUIRE( ic4.integrand(10) == 10+2-10-5);
    };

    SECTION ( " -= " ) {

        ic1 -= ic2;

        REQUIRE( ic1.number_of_integration_variables == 4);
        REQUIRE( ic1.integrand(10) == 10+2-10-5);

    };

};

TEST_CASE( "Binary operators + and - with functions which access private fields", "[IntegrandContainer]" ) {

    struct test_container1_t {
    private:
        const int test_field = 5;
    public:
        std::function<int(int)> func = [this] (int i) { return i+test_field; };
    } test_container1;

    struct test_container2_t {
    private:
        const int test_field = 9;
    public:
        std::function<int(int)> func = [this] (int i) { return i+test_field; };
    } test_container2;

    auto ic1 = secdecutil::IntegrandContainer<int, int>(3,test_container1.func);
    auto ic2 = secdecutil::IntegrandContainer<int, int>(4,test_container2.func);
    auto ic3 = ic1 + ic2;
    auto ic4 = ic1 - ic2;

    SECTION ( " + " ) {

        REQUIRE( ic3.number_of_integration_variables == 4);
        REQUIRE( ic3.integrand(10) == 10+5+10+9);

    };

    SECTION ( " += " ) {

        ic1 += ic2;

        REQUIRE( ic1.number_of_integration_variables == 4);
        REQUIRE( ic1.integrand(10) == 10+5+10+9);

    };

    SECTION ( " - " ) {

        REQUIRE( ic4.number_of_integration_variables == 4);
        REQUIRE( ic4.integrand(10) == 10+5-10-9);
    };

    SECTION ( " -= " ) {

        ic1 -= ic2;

        REQUIRE( ic1.number_of_integration_variables == 4);
        REQUIRE( ic1.integrand(10) == 10+5-10-9);

    };

};

TEST_CASE( "std::accumulate", "[IntegrandContainer]" ) {

    std::function<int(int)> func1 = [] (int i) { return i+2; };
    std::function<int(int)> func2 = [] (int i) { return i+5; };
    std::function<int(int)> func3 = [] (int i) { return i+7; };

    auto ic1 = secdecutil::IntegrandContainer<int, int>(3,func1);
    auto ic2 = secdecutil::IntegrandContainer<int, int>(4,func2);
    auto ic3 = secdecutil::IntegrandContainer<int, int>(5,func3);

    std::array<secdecutil::IntegrandContainer<int, int>, 3> ics = {ic1,ic2,ic3};

    SECTION ( "accumulate" ) {

        auto ic4 = std::accumulate(ics.begin()+1, ics.end(), ics.at(0));

        REQUIRE( ic4.number_of_integration_variables == 5);
        REQUIRE( ic4.integrand(10) == 10+2+10+5+10+7);

    };

}

TEST_CASE( "Unary operators + and -", "[IntegrandContainer]" ) {

    auto function = [] (int i, double d) { return i*d; };
    auto ic = secdecutil::IntegrandContainer<double, int, double>(2,function);

    SECTION ( " + " ) {

        REQUIRE( (+ic).number_of_integration_variables == 2 );
        REQUIRE( (+ic).integrand(3, 1.5) == Approx(4.5) );
        REQUIRE( (+ic).integrand(4, 1.5) == Approx(6.0) );

    };

    SECTION ( " - " ) {

        REQUIRE( (-ic).number_of_integration_variables == 2 );
        REQUIRE( (-ic).integrand(3, 1.5) == Approx(-4.5) );
        REQUIRE( (-ic).integrand(4, 1.5) == Approx(-6.0) );

    };

}
