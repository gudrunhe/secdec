#include "catch.hpp"
#include "../src/series.hpp"

#include <complex>
#include <sstream>

TEST_CASE( "Constructor exceptions for vector constructor", "[Series]" ) {

    // no content
    std::vector<int> vector_zero = {};
    REQUIRE_THROWS_AS( secdecutil::Series<int>(0,0,vector_zero), std::invalid_argument);

    // content size too small
    std::vector<int> vector_one = {1};
    REQUIRE_THROWS_AS( secdecutil::Series<int>(0,1,vector_one), std::invalid_argument);

    // content size too large
    std::vector<int> vector_three = {1,2,3};
    REQUIRE_THROWS_AS( secdecutil::Series<int>(0,1,vector_three), std::invalid_argument);

};

TEST_CASE( "Constructor exceptions for initializer list constructor", "[Series]" ) {

    // no content
    REQUIRE_THROWS_AS( secdecutil::Series<int>(0,0,{}), std::invalid_argument);

    // content size too small
    REQUIRE_THROWS_AS( secdecutil::Series<int>(0,1,{1}), std::invalid_argument);

    // content size too large
    REQUIRE_THROWS_AS( secdecutil::Series<int>(0,1,{1,2,3}), std::invalid_argument);

};

TEST_CASE( "Check Access", "[Series]" ) {

    int order_min = -2;
    int order_max = 1;
    bool truncated_above = false;

    // vector constructor
    std::vector<int> test_vector = {1,2,3,4};
    auto test_vector_original_value = test_vector.at(0);
    auto test_vector_series = secdecutil::Series<int>(order_min,order_max,test_vector,truncated_above);

    // initializer constructor
    auto test_list_series = secdecutil::Series<int>(order_min,order_max,{1,2,3,4},truncated_above);

    SECTION( "Accessing fields" ) {

        REQUIRE( test_vector_series.get_order_min() == order_min );
        REQUIRE( test_vector_series.get_order_max() == order_max );
        REQUIRE( test_vector_series.get_truncated_above() == truncated_above );

        REQUIRE( test_list_series.get_order_min() == order_min );
        REQUIRE( test_list_series.get_order_max() == order_max );
        REQUIRE( test_list_series.get_truncated_above() == truncated_above );

    };

    SECTION( "Accessing elements of series with [] and at()" ) {

        // operator[]
        REQUIRE( test_vector_series[-2] == test_vector.at(0) );
        REQUIRE( test_vector_series[-1] == test_vector.at(1) );
        REQUIRE( test_vector_series[0] == test_vector.at(2) );
        REQUIRE( test_vector_series[1] == test_vector.at(3) );

        REQUIRE( test_list_series[-2] == test_vector.at(0) );
        REQUIRE( test_list_series[-1] == test_vector.at(1) );
        REQUIRE( test_list_series[0] == test_vector.at(2) );
        REQUIRE( test_list_series[1] == test_vector.at(3) );

        // at()
        REQUIRE( test_vector_series.at(-2) == test_vector.at(0) );
        REQUIRE( test_vector_series.at(-1) == test_vector.at(1) );
        REQUIRE( test_vector_series.at(0) == test_vector.at(2) );
        REQUIRE( test_vector_series.at(1) == test_vector.at(3) );

        REQUIRE( test_list_series.at(-2) == test_vector.at(0) );
        REQUIRE( test_list_series.at(-1) == test_vector.at(1) );
        REQUIRE( test_list_series.at(0) == test_vector.at(2) );
        REQUIRE( test_list_series.at(1) == test_vector.at(3) );

    };

    SECTION( "Modifying elements with []" ) {

        test_vector_series[0] = -test_vector.at(2);

        test_list_series[0] = -test_vector.at(2);

        REQUIRE( test_vector_series.at(0) == -test_vector.at(2) );

        REQUIRE( test_list_series.at(0) == -test_vector.at(2) );

    };

    SECTION( "Modifying vector does not change series [] or at()" ) {

        // Change test_vector value
        test_vector.at(0) = test_vector.at(0) + 1.;

        // Series should have value contained in vector on construction, not currently
        REQUIRE( test_vector_series[-2] == test_vector_original_value );
        REQUIRE( test_vector_series.at(-2) == test_vector_original_value );

    };

    SECTION( "Modifying series does not change vector" ) {

        // Change series value
        test_vector_series[-2] = test_vector_series.at(-2) + 1.;

        // Vector should still have original value
        REQUIRE( test_vector.at(0) == test_vector_original_value );

    };

};

TEST_CASE( "Operators == and !=", "[Series]" ) {

    auto series_one = secdecutil::Series<int>(-3,4,{-5,-6,-7,0,2,3,4,-1});
    auto series_two = secdecutil::Series<int>(-3,4,{-5,-6,-7,0,2,3,4,-1});
    auto series_three = secdecutil::Series<int>(-3,4,{-5,-6,-7,0,2,3,4,-1},false);

    auto series_unidentical_one = secdecutil::Series<int>(-3,4,{-5,-6,-6,0,2,3,4,-1});
    auto series_unidentical_two = secdecutil::Series<int>(-4,4,{0,-5,-6,-7,0,2,3,4,-1});
    auto series_unidentical_three = secdecutil::Series<int>(-3,5,{-5,-6,-7,0,2,3,4,-1,2});

    REQUIRE( series_one == series_two );
    REQUIRE( series_one != series_unidentical_one );
    REQUIRE( series_one != series_unidentical_two );
    REQUIRE( series_one != series_unidentical_three );

};

TEST_CASE ( "Unary Operators - and +", "[Series]" ) {

    auto series_one = secdecutil::Series<int>(-3,4,{-5,-6,-7,0,2,3,4,-1});
    auto series_one_minus = secdecutil::Series<int>(-3,4,{5,6,7,0,-2,-3,-4,1});

    REQUIRE ( -series_one == series_one_minus );
    REQUIRE ( +series_one == series_one );

};

TEST_CASE( "Operator +", "[Series]" ) {

    auto series_one = secdecutil::Series<int>(                     -3,4,{-5 ,-6 ,-7   , 0   , 2   , 3  , 4   ,-1          });
    auto series_exact_one = secdecutil::Series<int>(               -3,4,{-5 ,-6 ,-7   , 0   , 2   , 3  , 4   ,-1          },false);
    auto series_two = secdecutil::Series<int>(                     -1,6,{         8   , 6   ,-7   , 0  , 8   , 5  , 4 , 1 });
    auto series_exact_two = secdecutil::Series<int>(               -1,6,{         8   , 6   ,-7   , 0  , 8   , 5  , 4 , 1 },false);
    auto series_one_plus_two = secdecutil::Series<int>(            -3,4,{-5 ,-6 ,-7+8 , 0+6 , 2-7 , 3+0, 4+8 ,-1+5        });
    auto series_one_plus_exact_two = secdecutil::Series<int>(      -3,4,{-5 ,-6 ,-7+8 , 0+6 , 2-7 , 3+0 ,4+8 ,-1+5        });
    auto series_exact_one_plus_two = secdecutil::Series<int>(      -3,6,{-5 ,-6 ,-7+8 , 0+6 , 2-7 , 3+0 ,4+8 ,-1+5, 4 , 1 });
    auto series_exact_one_plus_exact_two = secdecutil::Series<int>(-3,6,{-5 ,-6 ,-7+8 , 0+6 , 2-7 , 3+0, 4+8 ,-1+5, 4 , 1 },false);

    SECTION ( " + " ) {
        REQUIRE( ( series_one + series_two ) == series_one_plus_two );
        REQUIRE( ( series_one + series_exact_two ) == series_one_plus_exact_two );
        REQUIRE( ( series_exact_one + series_two ) == series_exact_one_plus_two );
        REQUIRE( ( series_exact_one + series_exact_two ) == series_exact_one_plus_exact_two );
    };

    SECTION ( "Test 1: +=" ) {
        REQUIRE( ( series_one += series_two ) == series_one_plus_two );
    };

    SECTION ( "Test 2: +=" ) {
        REQUIRE( ( series_one += series_exact_two ) == series_one_plus_exact_two );
    };

    SECTION ( "Test 3: +=" ) {
        REQUIRE( ( series_exact_one += series_two ) == series_exact_one_plus_two );
    };

    SECTION ( "Test 4: +=" ) {
        REQUIRE( ( series_exact_one += series_exact_two ) == series_exact_one_plus_exact_two );
    };

};

TEST_CASE( "Operator -", "[Series]" ) {

    auto series_exact_one = secdecutil::Series<int>(                       -3,4,{-5 ,-6 ,-7   , 0   , 2   , 3  , 4   ,-1          },false);
    auto series_exact_two = secdecutil::Series<int>(                       -1,6,{         8   , 6   ,-7   , 0  , 8   , 5  , 4 , 1 },false);
    auto series_exact_one_minus_series_exact_two = secdecutil::Series<int>(-3,6,{-5 ,-6 ,-7-8 , 0-6 , 2+7 , 3-0, 4-8 ,-1-5,-4 ,-1 },false);

    // Check behaviour for double as we use T(-1) in the - operator
    auto d_series_exact_one = secdecutil::Series<double>(                       -3,4,{-5. ,-6. ,-7.    , 0.    , 2.    , 3.   , 4.    ,-1.             },false);
    auto d_series_exact_two = secdecutil::Series<double>(                       -1,6,{           8.    , 6.    ,-7.    , 0.   , 8.    , 5.   , 4. , 1. },false);
    auto d_series_exact_one_minus_series_exact_two = secdecutil::Series<double>(-3,6,{-5. ,-6. ,-7.-8. , 0.-6. , 2.+7. , 3.-0., 4.-8. ,-1.-5.,-4. ,-1. },false);

    SECTION ( " - " ) {
        REQUIRE( ( series_exact_one - series_exact_two ) == series_exact_one_minus_series_exact_two );
        REQUIRE( ( d_series_exact_one - d_series_exact_two ) == d_series_exact_one_minus_series_exact_two );
    };

    SECTION ( "Test 1: -= " ) {
        REQUIRE( ( series_exact_one -= series_exact_two ) == series_exact_one_minus_series_exact_two );
    };

    SECTION ( "Test 2: -= " ) {
        REQUIRE( ( d_series_exact_one -= d_series_exact_two ) == d_series_exact_one_minus_series_exact_two );
    };

};

TEST_CASE( "Operator *", "[Series]" ) {

    auto series_one =  secdecutil::Series<int>(                 -2,1,{-5      ,-6      ,-7      , 3          });
    auto series_exact_one =  secdecutil::Series<int>(           -2,1,{-5      ,-6      ,-7      , 3          },false);
    auto series_two =  secdecutil::Series<int>(                 -1,3,{         -3      , 9      , 1   , 2, 3 });
    auto series_exact_two =  secdecutil::Series<int>(           -1,3,{         -3      , 9      , 1   , 2, 3 },false);
    auto three_times_series_exact_one = secdecutil::Series<int>(-2,1,{ 3*(-5) , 3*(-6) , 3*(-7) , 3*3        },false);
    auto series_exact_one_times_series_exact_two = secdecutil::Series<int>(-3,4,
                                                                           {
                                                                               (-5)*(-3),
                                                                               (-5)*9+(-6)*(-3),
                                                                               (-5)*1+(-6)*9+(-7)*(-3),
                                                                               (-5)*2+(-6)*1+(-7)*9+3*(-3),
                                                                               (-5)*3+(-6)*2+(-7)*1+3*9,
                                                                               (-6)*3+(-7)*2+3*1,
                                                                               (-7)*3+3*2,
                                                                               3*3
                                                                           },false);
    auto series_one_times_series_exact_two = secdecutil::Series<int>(-3,0,
                                                                     {
                                                                         (-5)*(-3),
                                                                         (-5)*9+(-6)*(-3),
                                                                         (-5)*1+(-6)*9+(-7)*(-3),
                                                                         (-5)*2+(-6)*1+(-7)*9+3*(-3)
                                                                     });
    auto series_exact_one_times_series_two = secdecutil::Series<int>(-3,1,
                                                                     {
                                                                         (-5)*(-3),
                                                                         (-5)*9+(-6)*(-3),
                                                                         (-5)*1+(-6)*9+(-7)*(-3),
                                                                         (-5)*2+(-6)*1+(-7)*9+3*(-3),
                                                                         (-5)*3+(-6)*2+(-7)*1+3*9
                                                                     });
    auto series_one_times_series_two = secdecutil::Series<int>(-3,0,
                                                               {
                                                                   (-5)*(-3),
                                                                   (-5)*9+(-6)*(-3),
                                                                   (-5)*1+(-6)*9+(-7)*(-3),
                                                                   (-5)*2+(-6)*1+(-7)*9+3*(-3)
                                                               });

    SECTION ( " * " ) {
        REQUIRE( ( series_exact_one * 3 ) == three_times_series_exact_one );
        REQUIRE( ( 3 * series_exact_one ) == three_times_series_exact_one );
        REQUIRE( ( series_exact_one * series_exact_two ) == series_exact_one_times_series_exact_two );
        REQUIRE( ( series_one * series_exact_two ) == series_one_times_series_exact_two );
        REQUIRE( ( series_exact_one * series_two ) == series_exact_one_times_series_two );
        REQUIRE( ( series_one * series_two ) == series_one_times_series_two );
    };

    SECTION ( "Test 1: *= " ) {
        REQUIRE( ( series_exact_one *= 3 ) == three_times_series_exact_one );
    };

    SECTION ( "Test 2: *= " ) {
        REQUIRE( ( series_exact_one *= series_exact_two ) == series_exact_one_times_series_exact_two );
    };

    SECTION ( "Test 3: *= " ) {
        REQUIRE( ( series_one *= series_exact_two ) == series_one_times_series_exact_two );
    };

    SECTION ( "Test 4: *= " ) {
        REQUIRE( ( series_exact_one *= series_two ) == series_exact_one_times_series_two );
    };

    SECTION ( "Test 5: *= " ) {
        REQUIRE( ( series_one *= series_two ) == series_one_times_series_two );
    };

};

TEST_CASE( "Operator * for complex", "[Series]" ) {

    auto one_plus_i_times_x = secdecutil::Series<std::complex<int>>(0,1,{{1,0},{0,1}});
    auto minus_one_minus_i_times_x = secdecutil::Series<std::complex<int>>(0,1,{{-1,0},{0,-1}});
    auto result_operator_minus = - one_plus_i_times_x;
    REQUIRE( result_operator_minus == minus_one_minus_i_times_x );

};

TEST_CASE( "Check Multivariate Access", "[Series]" ) {

    // general series
    auto multivariate_series =
    secdecutil::Series<secdecutil::Series<int>>(-1,1,{
        secdecutil::Series<int>(-1,0,{1,2}),
        secdecutil::Series<int>(0,2,{3,4,5}),
        secdecutil::Series<int>(2,5,{6,7,8,9})
    });

    SECTION( "Multivariate Accessing elements of series with [] and at()" ) {
        REQUIRE( multivariate_series[-1][-1] == 1 );
        REQUIRE( multivariate_series[-1][0] == 2 );
        REQUIRE( multivariate_series[0][0] == 3 );
        REQUIRE( multivariate_series[0][1] == 4 );
        REQUIRE( multivariate_series[0][2] == 5 );
        REQUIRE( multivariate_series[1][2] == 6 );
        REQUIRE( multivariate_series[1][3] == 7 );
        REQUIRE( multivariate_series[1][4] == 8 );
        REQUIRE( multivariate_series[1][5] == 9 );
    };

};

TEST_CASE( "Multivariate Operator +" ) {
    auto multivariate_series_one =
    secdecutil::Series<secdecutil::Series<int>>(0,2,{
        secdecutil::Series<int>(0,2,{1,1,1}),
        secdecutil::Series<int>(0,2,{1,1,1}),
        secdecutil::Series<int>(0,2,{1,1,1})
    });

    auto multivariate_series_two =
    secdecutil::Series<secdecutil::Series<int>>(-1,1,{
        secdecutil::Series<int>(-1,1,{1,1,1}),
        secdecutil::Series<int>(-1,1,{1,1,1}),
        secdecutil::Series<int>(-1,1,{1,1,1}),
    });

    // Note: result is not (-1,1,{1,2,2,1}) as last 1 is truncated
    auto multivariate_series_one_plus_two =
    secdecutil::Series<secdecutil::Series<int>>(-1,1,{
        secdecutil::Series<int>(-1,1,{1,1,1}),
        secdecutil::Series<int>(-1,1,{1,2,2}),
        secdecutil::Series<int>(-1,1,{1,2,2})
    });

    SECTION ( "+" ) {
        REQUIRE( (multivariate_series_one + multivariate_series_two) == multivariate_series_one_plus_two );
    };

    SECTION ( "Test 1: += " )
    {
        REQUIRE( (multivariate_series_one += multivariate_series_two) == multivariate_series_one_plus_two );
    }

};

TEST_CASE( "Multivariate Operator -" ) {
    auto multivariate_series_one =
    secdecutil::Series<secdecutil::Series<int>>(0,2,{
        secdecutil::Series<int>(0,2,{1,1,1}),
        secdecutil::Series<int>(0,2,{1,1,1}),
        secdecutil::Series<int>(0,2,{1,1,1})
    });

    auto multivariate_series_two =
    secdecutil::Series<secdecutil::Series<int>>(-1,1,{
        secdecutil::Series<int>(-1,1,{1,1,1}),
        secdecutil::Series<int>(-1,1,{1,1,1}),
        secdecutil::Series<int>(-1,1,{1,1,1}),
    });

    // Note: result is not (-1,1,{1,2,2,1}) as last 1 is truncated
    auto multivariate_series_one_minus_two =
    secdecutil::Series<secdecutil::Series<int>>(-1,1,{
        secdecutil::Series<int>(-1,1,{-1,-1,-1}),
        secdecutil::Series<int>(-1,1,{-1,0,0}),
        secdecutil::Series<int>(-1,1,{-1,0,0})
    });

    SECTION ( "-" ) {
        REQUIRE( (multivariate_series_one - multivariate_series_two) == multivariate_series_one_minus_two );
    };

    SECTION ( "Test 1: -= " )
    {
        REQUIRE( (multivariate_series_one -= multivariate_series_two) == multivariate_series_one_minus_two );
    }

};

TEST_CASE( "Multivariate Operator *" ) {

    auto multivariate_series_one =
    secdecutil::Series<secdecutil::Series<int>>(0,1,{
        secdecutil::Series<int>(-1,0,{1,1},false),
        secdecutil::Series<int>(0,1,{1,1},false),
    },false);

    auto multivariate_series_one_sq =
    secdecutil::Series<secdecutil::Series<int>>(0,2,{
        secdecutil::Series<int>(-2,0,{1,2,1},false),
        secdecutil::Series<int>(-1,1,{2,4,2},false),
        secdecutil::Series<int>(0,2,{1,2,1},false)
    },false);

    SECTION ( "*" ) {
        REQUIRE( (multivariate_series_one * multivariate_series_one) == multivariate_series_one_sq );
    };

    SECTION ( "Test 1: *=" ) {
        REQUIRE( (multivariate_series_one *= multivariate_series_one) == multivariate_series_one_sq );
    };
};

TEST_CASE( "Naming the expansion parameter for operator <<" ) {

    // default expansion parameter: 'x'
    auto one_plus_x = secdecutil::Series<int>(0,1,{1,1},false);
    std::stringstream stream_x;
    stream_x << one_plus_x;
    std::string target_x = " + (1) + (1)*x";
    REQUIRE( stream_x.str() == target_x );

    auto one_plus_eps = secdecutil::Series<int>(0,2,{1,1,8},true);
    one_plus_eps.expansion_parameter = "eps"; // rename
    std::stringstream stream_eps;
    stream_eps << one_plus_eps;
    std::string target_eps = " + (1) + (1)*eps + (8)*eps^2 + O(eps^3)";
    REQUIRE( stream_eps.str() == target_eps );

    auto plus_order_one = secdecutil::Series<int>(-3,-1,{1,1,8},true);
    std::stringstream stream_order_one;
    stream_order_one << plus_order_one;
    std::string target_order_one = " + (1)*x^-3 + (1)*x^-2 + (8)*x^-1 + O(1)";
    REQUIRE( stream_order_one.str() == target_order_one );

    auto plus_order_eps = secdecutil::Series<int>(-2,0,{1,1,8},true);
    plus_order_eps.expansion_parameter = "eps"; // rename
    std::stringstream stream_order_eps;
    stream_order_eps << plus_order_eps;
    std::string target_order_eps = " + (1)*eps^-2 + (1)*eps^-1 + (8) + O(eps)";
    REQUIRE( stream_order_eps.str() == target_order_eps );

};
