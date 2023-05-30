#include <vector>
#include <functional>

#include "../secdecutil/deep_apply.hpp"
#include "../secdecutil/series.hpp"

#include "catch_amalgamated.hpp"
using Catch::Approx;

TEST_CASE( "deep_apply std::vector container to different type", "[deep_apply]" ) {

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

        auto zero_nest_double_to_int = secdecutil::deep_apply(zero_nest_double, double_to_int);

        REQUIRE( typeid(zero_nest_double_to_int) == typeid(int) );
        REQUIRE( zero_nest_double_to_int == zero_nest_int );

    };

    SECTION( "type conversion 1-nest" ) {

        auto one_nest_double_to_int = secdecutil::deep_apply(one_nest_double, double_to_int);

        REQUIRE( typeid(one_nest_double_to_int) == typeid(one_nest_int) );
        REQUIRE( one_nest_double_to_int == one_nest_int );

    };

    SECTION( "type conversion 2-nest" ) {

        auto two_nest_double_to_int = secdecutil::deep_apply(two_nest_double, double_to_int);

        REQUIRE( typeid(two_nest_double_to_int) == typeid(two_nest_int) );
        REQUIRE( two_nest_double_to_int == two_nest_int );

    };

    SECTION( "type conversion 3-nest" ) {

        auto three_nest_double_to_int = secdecutil::deep_apply(three_nest_double, double_to_int);

        REQUIRE( typeid(three_nest_double_to_int) == typeid(three_nest_int) );
        REQUIRE( three_nest_double_to_int == three_nest_int );

    };

};

TEST_CASE( "deep_apply const std::vector container to different type", "[deep_apply]" ) {

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

        auto zero_nest_double_to_int = secdecutil::deep_apply(zero_nest_double, double_to_int);

        REQUIRE( typeid(zero_nest_double_to_int) == typeid(int) );
        REQUIRE( zero_nest_double_to_int == zero_nest_int );

    };

    SECTION( "type conversion 1-nest" ) {

        auto one_nest_double_to_int = secdecutil::deep_apply(one_nest_double, double_to_int);

        REQUIRE( typeid(one_nest_double_to_int) == typeid(one_nest_int) );
        REQUIRE( one_nest_double_to_int == one_nest_int );

    };

    SECTION( "type conversion 2-nest" ) {

        auto two_nest_double_to_int = secdecutil::deep_apply(two_nest_double, double_to_int);

        REQUIRE( typeid(two_nest_double_to_int) == typeid(two_nest_int) );
        REQUIRE( two_nest_double_to_int == two_nest_int );

    };


    SECTION( "type conversion 3-nest" ) {

        auto three_nest_double_to_int = secdecutil::deep_apply(three_nest_double, double_to_int);

        REQUIRE( typeid(three_nest_double_to_int) == typeid(three_nest_int) );
        REQUIRE( three_nest_double_to_int == three_nest_int );

    };

};

TEST_CASE( "increment std::vector container", "[deep_apply]" ) {

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
    const std::function<void(int&)> increment = [] (int& a) { a++; };

    SECTION( "increment 0-nest and return 0" ) {

        auto zero_nest_int_to_zero = secdecutil::deep_apply(zero_nest_int, increment_and_return_zero);

        REQUIRE( zero_nest_int == zero_nest_int_increment );
        REQUIRE( zero_nest_int_to_zero == zero_nest_int_zero );

    };

    SECTION( "increment 1-nest and return 0" ) {

        auto one_nest_int_to_zero = secdecutil::deep_apply(one_nest_int, increment_and_return_zero);

        REQUIRE( one_nest_int == one_nest_int_increment );
        REQUIRE( one_nest_int_to_zero == one_nest_int_zero );

    };

    SECTION( "increment 2-nest and return 0" ) {

        auto two_nest_int_to_zero = secdecutil::deep_apply(two_nest_int, increment_and_return_zero);

        REQUIRE( two_nest_int == two_nest_int_increment );
        REQUIRE( two_nest_int_to_zero == two_nest_int_zero );

    };

    SECTION( "increment 3-nest and return 0" ) {

        auto three_nest_int_to_zero = secdecutil::deep_apply(three_nest_int, increment_and_return_zero);

        REQUIRE( three_nest_int == three_nest_int_increment );
        REQUIRE( three_nest_int_to_zero == three_nest_int_zero );

    };

    SECTION( "increment 0-nest" ) {

        secdecutil::deep_apply(zero_nest_int, increment);

        REQUIRE( zero_nest_int == zero_nest_int_increment );

    };

    SECTION( "increment 1-nest" ) {

        secdecutil::deep_apply(one_nest_int, increment);

        REQUIRE( one_nest_int == one_nest_int_increment );

    };

    SECTION( "increment 2-nest" ) {

        secdecutil::deep_apply(two_nest_int, increment);

        REQUIRE( two_nest_int == two_nest_int_increment );

    };

    SECTION( "increment 3-nest" ) {

        secdecutil::deep_apply(three_nest_int, increment);

        REQUIRE( three_nest_int == three_nest_int_increment );

    };

};

TEST_CASE( "deep_apply secdecutil::Series container to different type", "[deep_apply]" ) {

    secdecutil::Series<double> one_nest_double = {-1,1,{1.0, 2.0, 3.0}};
    secdecutil::Series<int>  one_nest_int = {-1,1,{1, 2, 3}};
    secdecutil::Series<secdecutil::Series<double>> two_nest_double =
    {0,1,
        {
            {-1,1, { 1.0, 2.0, 3.0}},
            {-2,0, { 4.0, 5.0, 6.0}}
        }
    };
    secdecutil::Series<secdecutil::Series<int>> two_nest_int =
    {0,1,
        {
            {-1,1, { 1, 2, 3}},
            {-2,0, { 4, 5, 6}}
        }
    };
    secdecutil::Series<secdecutil::Series<secdecutil::Series<double>>> three_nest_double =
    {0,1,
        {
            {-1,1,
                {
                    {-1,1, { 1.0, 2.0, 3.0}},
                    {-2,1, { 4.0, 5.0, 6.0, 7.0}},
                    { 1,3, { 8.0, 9.0, 10.0}}
                }
            },
            {0,0,
                {
                    {-1,2, { 11.0, 12.0, 13.0, 14.0}}
                }
            }
        }
    };
    secdecutil::Series<secdecutil::Series<secdecutil::Series<int>>> three_nest_int =
    {0,1,
        {
            {-1,1,
                {
                    {-1,1, { 1, 2, 3}},
                    {-2,1, { 4, 5, 6, 7}},
                    { 1,3, { 8, 9, 10}}
                }
            },
            {0,0,
                {
                    {-1,2, { 11, 12, 13, 14}}
                }
            }
        }
    };

    const std::function<int(double)> double_to_int = [] (double a) { return static_cast<int>(a); };

    SECTION( "type conversion 1-nest" ) {

        auto one_nest_double_to_int = secdecutil::deep_apply(one_nest_double, double_to_int);

        REQUIRE( typeid(one_nest_double_to_int) == typeid(one_nest_int) );
        REQUIRE( one_nest_double_to_int == one_nest_int );

    };

    SECTION( "type conversion 2-nest" ) {

        auto two_nest_double_to_int = secdecutil::deep_apply(two_nest_double, double_to_int);

        REQUIRE( typeid(two_nest_double_to_int) == typeid(two_nest_int) );
        REQUIRE( two_nest_double_to_int == two_nest_int );

    };


    SECTION( "type conversion 3-nest" ) {

        auto three_nest_double_to_int = secdecutil::deep_apply(three_nest_double, double_to_int);

        REQUIRE( typeid(three_nest_double_to_int) == typeid(three_nest_int) );
        REQUIRE( three_nest_double_to_int == three_nest_int );

    };

};

TEST_CASE( "deep_apply const secdecutil::Series container to different type", "[deep_apply]" ) {

    const secdecutil::Series<double> one_nest_double = {-1,1,{1.0, 2.0, 3.0}};
    const secdecutil::Series<int>  one_nest_int = {-1,1,{1, 2, 3}};
    const secdecutil::Series<secdecutil::Series<double>> two_nest_double =
    {0,1,
        {
            {-1,1, { 1.0, 2.0, 3.0}},
            {-2,0, { 4.0, 5.0, 6.0}}
        }
    };
    const secdecutil::Series<secdecutil::Series<int>> two_nest_int =
    {0,1,
        {
            {-1,1, { 1, 2, 3}},
            {-2,0, { 4, 5, 6}}
        }
    };
    const secdecutil::Series<secdecutil::Series<secdecutil::Series<double>>> three_nest_double =
    {0,1,
        {
            {-1,1,
                {
                    {-1,1, { 1.0, 2.0, 3.0}},
                    {-2,1, { 4.0, 5.0, 6.0, 7.0}},
                    { 1,3, { 8.0, 9.0, 10.0}}
                }
            },
            {0,0,
                {
                    {-1,2, { 11.0, 12.0, 13.0, 14.0}}
                }
            }
        }
    };
    const secdecutil::Series<secdecutil::Series<secdecutil::Series<int>>> three_nest_int =
    {0,1,
        {
            {-1,1,
                {
                    {-1,1, { 1, 2, 3}},
                    {-2,1, { 4, 5, 6, 7}},
                    { 1,3, { 8, 9, 10}}
                }
            },
            {0,0,
                {
                    {-1,2, { 11, 12, 13, 14}}
                }
            }
        }
    };

    const std::function<int(double)> double_to_int = [] (double a) { return static_cast<int>(a); };

    SECTION( "type conversion 1-nest" ) {

        auto one_nest_double_to_int = secdecutil::deep_apply(one_nest_double, double_to_int);

        REQUIRE( typeid(one_nest_double_to_int) == typeid(one_nest_int) );
        REQUIRE( one_nest_double_to_int == one_nest_int );

    };

    SECTION( "type conversion 2-nest" ) {

        auto two_nest_double_to_int = secdecutil::deep_apply(two_nest_double, double_to_int);

        REQUIRE( typeid(two_nest_double_to_int) == typeid(two_nest_int) );
        REQUIRE( two_nest_double_to_int == two_nest_int );

    };


    SECTION( "type conversion 3-nest" ) {

        auto three_nest_double_to_int = secdecutil::deep_apply(three_nest_double, double_to_int);

        REQUIRE( typeid(three_nest_double_to_int) == typeid(three_nest_int) );
        REQUIRE( three_nest_double_to_int == three_nest_int );

    };

};

TEST_CASE( "increment secdecutil::Series container", "[deep_apply]" ) {


    secdecutil::Series<int>  one_nest_int = {-1,1,{1, 2, 3}};
    secdecutil::Series<int>  one_nest_int_increment = {-1,1,{2, 3, 4}};
    secdecutil::Series<int>  one_nest_int_zero = {-1,1,{0, 0, 0}};
    secdecutil::Series<secdecutil::Series<int>> two_nest_int =
    {0,1,
        {
            {-1,1, { 1, 2, 3}},
            {-2,0, { 4, 5, 6}}
        }
    };
    secdecutil::Series<secdecutil::Series<int>> two_nest_int_increment =
    {0,1,
        {
            {-1,1, { 2, 3, 4}},
            {-2,0, { 5, 6, 7}}
        }
    };
    secdecutil::Series<secdecutil::Series<int>> two_nest_int_zero =
    {0,1,
        {
            {-1,1, { 0, 0, 0}},
            {-2,0, { 0, 0, 0}}
        }
    };
    secdecutil::Series<secdecutil::Series<secdecutil::Series<int>>> three_nest_int =
    {0,1,
        {
            {-1,1,
                {
                    {-1,1, { 1, 2, 3}},
                    {-2,1, { 4, 5, 6, 7}},
                    { 1,3, { 8, 9, 10}}
                }
            },
            {0,0,
                {
                    {-1,2, { 11, 12, 13, 14}}
                }
            }
        }
    };
    secdecutil::Series<secdecutil::Series<secdecutil::Series<int>>> three_nest_int_increment =
    {0,1,
        {
            {-1,1,
                {
                    {-1,1, { 2, 3, 4}},
                    {-2,1, { 5, 6, 7, 8}},
                    { 1,3, { 9, 10, 11}}
                }
            },
            {0,0,
                {
                    {-1,2, { 12, 13, 14, 15}}
                }
            }
        }
    };
    secdecutil::Series<secdecutil::Series<secdecutil::Series<int>>> three_nest_int_zero =
    {0,1,
        {
            {-1,1,
                {
                    {-1,1, { 0, 0, 0}},
                    {-2,1, { 0, 0, 0, 0}},
                    { 1,3, { 0, 0, 0}}
                }
            },
            {0,0,
                {
                    {-1,2, { 0, 0, 0, 0}}
                }
            }
        }
    };

    const std::function<int(int&)> increment_and_return_zero = [] (int& a) { a++; return 0; };
    const std::function<void(int&)> increment = [] (int& a) { a++; };

    SECTION( "increment 1-nest and return 0" ) {

        auto one_nest_int_to_zero = secdecutil::deep_apply(one_nest_int, increment_and_return_zero);

        REQUIRE( one_nest_int == one_nest_int_increment );
        REQUIRE( one_nest_int_to_zero == one_nest_int_zero );

    };

    SECTION( "increment 2-nest and return 0" ) {

        auto two_nest_int_to_zero = secdecutil::deep_apply(two_nest_int, increment_and_return_zero);

        REQUIRE( two_nest_int == two_nest_int_increment );
        REQUIRE( two_nest_int_to_zero == two_nest_int_zero );

    };

    SECTION( "increment 3-nest and return 0" ) {

        auto three_nest_int_to_zero = secdecutil::deep_apply(three_nest_int, increment_and_return_zero);

        REQUIRE( three_nest_int == three_nest_int_increment );
        REQUIRE( three_nest_int_to_zero == three_nest_int_zero );

    };

    SECTION( "increment 1-nest" ) {

        secdecutil::deep_apply(one_nest_int, increment);

        REQUIRE( one_nest_int == one_nest_int_increment );

    };

    SECTION( "increment 2-nest" ) {

        secdecutil::deep_apply(two_nest_int, increment);

        REQUIRE( two_nest_int == two_nest_int_increment );

    };

    SECTION( "increment 3-nest" ) {

        secdecutil::deep_apply(three_nest_int, increment);

        REQUIRE( three_nest_int == three_nest_int_increment );

    };

};

TEST_CASE( "side-effect only function on const std::vector container", "[deep_apply]" ) {

    const int zero_nest_int = 1;
    const std::vector<int> one_nest_int = { 1, 2, 3};

    const std::vector<std::vector<int>> two_nest_int =
    {
        { 1, 2, 3},
        { 4, 5, 6}
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

    const std::function<void(const int&, int&)> add_one_to_out = [] (const int& in, int& out) { ++out; };

    SECTION( "count 0-nest" ) {

        int number_of_elements = 0;
        const std::function<void(const int&)> func = std::bind(add_one_to_out, std::placeholders::_1, std::ref(number_of_elements));
        secdecutil::deep_apply(zero_nest_int, func);

        REQUIRE( number_of_elements == 1 );

    };

    SECTION( "count 1-nest" ) {

        int number_of_elements = 0;
        const std::function<void(const int&)> func = std::bind(add_one_to_out, std::placeholders::_1, std::ref(number_of_elements));
        secdecutil::deep_apply(one_nest_int, func);

        REQUIRE( number_of_elements == 3 );
    };

    SECTION( "count 2-nest" ) {

        int number_of_elements = 0;
        const std::function<void(const int&)> func = std::bind(add_one_to_out, std::placeholders::_1, std::ref(number_of_elements));
        secdecutil::deep_apply(two_nest_int, func);

        REQUIRE( number_of_elements == 6 );
    };

    SECTION( "count 3-nest" ) {

        int number_of_elements = 0;
        const std::function<void(const int&)> func = std::bind(add_one_to_out, std::placeholders::_1, std::ref(number_of_elements));
        secdecutil::deep_apply(three_nest_int, func);

        REQUIRE( number_of_elements == 14 );

    };

};

TEST_CASE( "act-on-innermost-container function", "[deep_apply]" ) {

    std::vector<std::vector<int>> two_nest_vector =
    {
            { 1, 2, 3},
            { 4, 5, 6, 7},
            { 8, 9, 10},
            { 11, 12, 13, 14}
    };

    secdecutil::Series<secdecutil::Series<int>> two_nest_series =
    {0,1,
        {
            {-1, 1, { 1, 1, 1}},
            {-2,-1, { 1, 1   },false /* not truncated */}
        }
    };

    SECTION("const void (side-effect only)") {

        size_t innermost_size;
        const std::function<void(const std::vector<int>&)> add_to_innermost_size_vector = [ &innermost_size ] (const std::vector<int>& v) { innermost_size += v.size(); };
        const std::function<void(const secdecutil::Series<int>&)> add_to_innermost_size_series = [ &innermost_size ] (const secdecutil::Series<int>& s) { innermost_size += s.get_content().size(); };

        innermost_size = 0;
        secdecutil::deep_apply(two_nest_vector, add_to_innermost_size_vector);
        REQUIRE( innermost_size == 14 );

        innermost_size = 0;
        secdecutil::deep_apply(two_nest_series, add_to_innermost_size_series);
        REQUIRE( innermost_size == 5 );

    };

    SECTION("non-const void (modify)") {

        const std::function<void(std::vector<int>&)> append_42_to_vector = [  ] (std::vector<int>& v) { v.push_back(42); };
        const std::function<void(secdecutil::Series<int>&)> add_42_to_series = [  ] (secdecutil::Series<int>& s) { s += 42; };

        secdecutil::deep_apply(two_nest_vector, append_42_to_vector);
        std::vector<std::vector<int>> target_two_nest_vector =
        {
            { 1, 2, 3, 42},
            { 4, 5, 6, 7, 42},
            { 8, 9, 10, 42},
            { 11, 12, 13, 14, 42}
        };
        REQUIRE(two_nest_vector.size() == 4);
        for(size_t i = 0; i < 4; ++i)
            REQUIRE(two_nest_vector.at(i) == target_two_nest_vector.at(i));

        secdecutil::deep_apply(two_nest_series, add_42_to_series);
        secdecutil::Series<secdecutil::Series<int>> target_two_nest_series =
        {0,1,
            {
                {-1,1, { 1, 43, 1}},
                {-2,0, { 1, 1, 42},false /* not truncated */}
            }
        };
        REQUIRE(two_nest_series == target_two_nest_series);

    };

    SECTION("const non-void (return)") {

        const std::function<size_t(const std::vector<int>&)> get_innermost_sizes_vector = [  ] (const std::vector<int>& v) { return v.size(); };
        const std::function<size_t(const secdecutil::Series<int>&)> get_innermost_sizes_series = [  ] (const secdecutil::Series<int>& s) { return s.get_content().size(); };

        auto sizes_two_nest_vector = secdecutil::deep_apply(two_nest_vector, get_innermost_sizes_vector);
        std::vector<size_t> target_sizes_two_nest_vector = {3,4,3,4};
        REQUIRE(typeid(sizes_two_nest_vector) == typeid(target_sizes_two_nest_vector));
        REQUIRE(sizes_two_nest_vector == target_sizes_two_nest_vector);

        auto sizes_two_nest_series = secdecutil::deep_apply(two_nest_series, get_innermost_sizes_series);
        secdecutil::Series<size_t> target_sizes_two_nest_series = {0,1,{3,2}};
        REQUIRE(typeid(sizes_two_nest_series) == typeid(target_sizes_two_nest_series));
        REQUIRE(sizes_two_nest_series == target_sizes_two_nest_series);

    };

    SECTION("non-const non-void (modify and return size_count)") {

        const std::function<size_t(std::vector<int>&)> append_42_to_vector = [  ] (std::vector<int>& v) { v.push_back(42); return v.size(); };
        const std::function<size_t(secdecutil::Series<int>&)> add_42_to_series = [  ] (secdecutil::Series<int>& s) { s += 42; return s.get_content().size(); };

        auto sizes_two_nest_vector = secdecutil::deep_apply(two_nest_vector, append_42_to_vector);
        std::vector<std::vector<int>> target_two_nest_vector =
        {
            { 1, 2, 3, 42},
            { 4, 5, 6, 7, 42},
            { 8, 9, 10, 42},
            { 11, 12, 13, 14, 42}
        };
        REQUIRE(typeid(sizes_two_nest_vector) == typeid(std::vector<size_t>));
        REQUIRE(two_nest_vector.size() == 4);
        REQUIRE(sizes_two_nest_vector.size() == 4);
        for(size_t i = 0; i < 4; ++i)
        {
            REQUIRE(two_nest_vector.at(i) == target_two_nest_vector.at(i));
            REQUIRE(sizes_two_nest_vector.at(i) == target_two_nest_vector.at(i).size());
        };

        auto sizes_two_nest_series = secdecutil::deep_apply(two_nest_series, add_42_to_series);
        secdecutil::Series<secdecutil::Series<int>> target_two_nest_series =
        {0,1,
            {
                {-1,1, { 1, 43, 1}},
                {-2,0, { 1, 1, 42},false /* not truncated */}
            }
        };
        secdecutil::Series<size_t> target_sizes_two_nest_series = {0,1,{3,3}};
        REQUIRE(two_nest_series == target_two_nest_series);
        REQUIRE(typeid(sizes_two_nest_series) == typeid(target_sizes_two_nest_series));
        REQUIRE(sizes_two_nest_series == target_sizes_two_nest_series);

    };

};
