#include <vector>
#include <functional>
#include "../secdecutil/algorithm.hpp"

#include "catch.hpp"

TEST_CASE( "transform const std::vector container to different type", "[transform]" ) {

    const double zero_nest_double = 1.0;
    const int zero_nest_int = 1;

    const std::vector<double> one_nest_double = { 1.0, 2.0, 3.0 };
    const std::vector<int> one_nest_int = { 1, 2, 3};

    const std::vector<std::vector<double>> two_nest_double =
    {
        { 1.0, 2.0, 3.0},
        { 4.0, 5.0, 6.0}
    };
    const std::vector<std::vector<int>> two_nest_int =
    {
        { 1, 2, 3},
        { 4, 5, 6}
    };

    const std::vector<std::vector<std::vector<double>>> three_nest_double =
    {
        {
            { 1.0, 2.0, 3.0},
            { 4.0, 5.0, 6.0, 7.0},
            { 8.0, 9.0, 10.0}
        },
        {
            { 11.0, 12.0, 13.0, 14.0}
        }
    };
    const std::vector<std::vector<std::vector<int>>> three_nest_int =
    {
        {
            { 1, 2, 3},
            { 4, 5, 6, 7},
            { 8, 9, 10}
        },
        {
            { 11, 12, 13, 14}
        }
    };

    const std::function<int(double)> double_to_int = [] (double a) { return static_cast<int>(a); };

    SECTION( "type conversion 0-nest" ) {

        auto zero_nest_double_to_int = secdecutil::transform(zero_nest_double, double_to_int);

        REQUIRE( typeid(zero_nest_double_to_int) == typeid(int) );
        REQUIRE( zero_nest_double_to_int == zero_nest_int );

    };

    SECTION( "type conversion 1-nest" ) {

        auto one_nest_double_to_int = secdecutil::transform(one_nest_double, double_to_int);

        REQUIRE( typeid(one_nest_double_to_int) == typeid(one_nest_int) );
        REQUIRE( one_nest_double_to_int == one_nest_int );

    };

    SECTION( "type conversion 2-nest" ) {

        auto two_nest_double_to_int = secdecutil::transform(two_nest_double, double_to_int);

        REQUIRE( typeid(two_nest_double_to_int) == typeid(two_nest_int) );
        REQUIRE( two_nest_double_to_int == two_nest_int );

    };


    SECTION( "type conversion 3-nest" ) {

        auto three_nest_double_to_int = secdecutil::transform(three_nest_double, double_to_int);

        REQUIRE( typeid(three_nest_double_to_int) == typeid(three_nest_int) );
        REQUIRE( three_nest_double_to_int == three_nest_int );

    };

};

TEST_CASE( "transform std::vector container to different type", "[transform]" ) {

    double zero_nest_double = 1.0;
    int zero_nest_int = 1;

    std::vector<double> one_nest_double = { 1.0, 2.0, 3.0 };
    std::vector<int> one_nest_int = { 1, 2, 3};

    std::vector<std::vector<double>> two_nest_double =
    {
        { 1.0, 2.0, 3.0},
        { 4.0, 5.0, 6.0}
    };
    std::vector<std::vector<int>> two_nest_int =
    {
        { 1, 2, 3},
        { 4, 5, 6}
    };

    std::vector<std::vector<std::vector<double>>> three_nest_double =
    {
        {
            { 1.0, 2.0, 3.0},
            { 4.0, 5.0, 6.0, 7.0},
            { 8.0, 9.0, 10.0}
        },
        {
            { 11.0, 12.0, 13.0, 14.0}
        }
    };
    std::vector<std::vector<std::vector<int>>> three_nest_int =
    {
        {
            { 1, 2, 3},
            { 4, 5, 6, 7},
            { 8, 9, 10}
        },
        {
            { 11, 12, 13, 14}
        }
    };

    const std::function<int(double)> double_to_int = [] (double a) { return static_cast<int>(a); };

    SECTION( "type conversion 0-nest" ) {

        auto zero_nest_double_to_int = secdecutil::transform(zero_nest_double, double_to_int);

        REQUIRE( typeid(zero_nest_double_to_int) == typeid(int) );
        REQUIRE( zero_nest_double_to_int == zero_nest_int );

    };

    SECTION( "type conversion 1-nest" ) {

        auto one_nest_double_to_int = secdecutil::transform(one_nest_double, double_to_int);

        REQUIRE( typeid(one_nest_double_to_int) == typeid(one_nest_int) );
        REQUIRE( one_nest_double_to_int == one_nest_int );

    };

    SECTION( "type conversion 2-nest" ) {

        auto two_nest_double_to_int = secdecutil::transform(two_nest_double, double_to_int);

        REQUIRE( typeid(two_nest_double_to_int) == typeid(two_nest_int) );
        REQUIRE( two_nest_double_to_int == two_nest_int );

    };


    SECTION( "type conversion 3-nest" ) {

        auto three_nest_double_to_int = secdecutil::transform(three_nest_double, double_to_int);

        REQUIRE( typeid(three_nest_double_to_int) == typeid(three_nest_int) );
        REQUIRE( three_nest_double_to_int == three_nest_int );
        
    };
    
};

TEST_CASE( "increment std::vector container and return 0", "[transform]" ) {


    int zero_nest_int = 1;
    int zero_nest_int_increment = 2;
    int zero_nest_int_zero = 0;

    std::vector<int> one_nest_int = { 1, 2, 3};
    std::vector<int> one_nest_int_increment = { 2, 3, 4};
    std::vector<int> one_nest_int_zero = { 0, 0, 0};

    std::vector<std::vector<int>> two_nest_int =
    {
        { 1, 2, 3},
        { 4, 5, 6}
    };
    std::vector<std::vector<int>> two_nest_int_increment =
    {
        { 2, 3, 4},
        { 5, 6, 7}
    };
    std::vector<std::vector<int>> two_nest_int_zero =
    {
        { 0, 0, 0},
        { 0, 0, 0}
    };

    std::vector<std::vector<std::vector<int>>> three_nest_int =
    {
        {
            { 1, 2, 3},
            { 4, 5, 6, 7},
            { 8, 9, 10}
        },
        {
            { 11, 12, 13, 14}
        }
    };
    std::vector<std::vector<std::vector<int>>> three_nest_int_increment =
    {
        {
            { 2, 3, 4},
            { 5, 6, 7, 8},
            { 9, 10, 11}
        },
        {
            { 12, 13, 14, 15}
        }
    };
    std::vector<std::vector<std::vector<int>>> three_nest_int_zero =
    {
        {
            { 0, 0, 0},
            { 0, 0, 0, 0},
            { 0, 0, 0}
        },
        {
            { 0, 0, 0, 0}
        }
    };

    const std::function<int(int&)> increment_and_return_zero = [] (int& a) { a++; return 0; };

    SECTION( "increment 0-nest" ) {

        auto zero_nest_int_to_zero = secdecutil::transform(zero_nest_int, increment_and_return_zero);

        REQUIRE( zero_nest_int == zero_nest_int_increment );
        REQUIRE( zero_nest_int_to_zero == zero_nest_int_zero );

    };

    SECTION( "increment 1-nest" ) {

        auto one_nest_int_to_zero = secdecutil::transform(one_nest_int, increment_and_return_zero);

        REQUIRE( one_nest_int == one_nest_int_increment );
        REQUIRE( one_nest_int_to_zero == one_nest_int_zero );

    };

    SECTION( "increment 2-nest" ) {

        auto two_nest_int_to_zero = secdecutil::transform(two_nest_int, increment_and_return_zero);

        REQUIRE( two_nest_int == two_nest_int_increment );
        REQUIRE( two_nest_int_to_zero == two_nest_int_zero );

    };

    SECTION( "increment 3-nest" ) {

        auto three_nest_int_to_zero = secdecutil::transform(three_nest_int, increment_and_return_zero);

        REQUIRE( three_nest_int == three_nest_int_increment );
        REQUIRE( three_nest_int_to_zero == three_nest_int_zero );

    };

};

TEST_CASE( "increment std::vector container", "[transform]" ) {


    int zero_nest_int = 1;
    int zero_nest_int_increment = 2;

    std::vector<int> one_nest_int = { 1, 2, 3};
    std::vector<int> one_nest_int_increment = { 2, 3, 4};

    std::vector<std::vector<int>> two_nest_int =
    {
        { 1, 2, 3},
        { 4, 5, 6}
    };
    std::vector<std::vector<int>> two_nest_int_increment =
    {
        { 2, 3, 4},
        { 5, 6, 7}
    };

    std::vector<std::vector<std::vector<int>>> three_nest_int =
    {
        {
            { 1, 2, 3},
            { 4, 5, 6, 7},
            { 8, 9, 10}
        },
        {
            { 11, 12, 13, 14}
        }
    };
    std::vector<std::vector<std::vector<int>>> three_nest_int_increment =
    {
        {
            { 2, 3, 4},
            { 5, 6, 7, 8},
            { 9, 10, 11}
        },
        {
            { 12, 13, 14, 15}
        }
    };

    const std::function<void(int&)> increment = [] (int& a) { a++; };

    SECTION( "increment 0-nest" ) {

        secdecutil::transform(zero_nest_int, increment);

        REQUIRE( zero_nest_int == zero_nest_int_increment );

    };

    SECTION( "increment 1-nest" ) {

        secdecutil::transform(one_nest_int, increment);

        REQUIRE( one_nest_int == one_nest_int_increment );

    };

    SECTION( "increment 2-nest" ) {

        secdecutil::transform(two_nest_int, increment);

        REQUIRE( two_nest_int == two_nest_int_increment );

    };

    SECTION( "increment 3-nest" ) {

        secdecutil::transform(three_nest_int, increment);

        REQUIRE( three_nest_int == three_nest_int_increment );

    };
    
};
