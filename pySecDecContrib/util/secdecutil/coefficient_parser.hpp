#include <cassert> // assert
#include <cmath> // std::frexp
#include <istream> // std::istream
#include <sstream> // std::stringstream
#include <string> // std::string
#include <vector> // std::vector

#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include <exparse.hpp>

namespace secdecutil
{
    namespace exparse
    {
    
        struct unknown_expression_error : public std::runtime_error { using std::runtime_error::runtime_error; };
    
        //namespace g = GiNaC;
        typedef std::vector<long long int> powerlist_t;
        typedef mpqc_class rational_t;
        typedef std::map<powerlist_t, rational_t> expression_t;

        const static std::string blankchars(" \t\n\v\f\r");
        inline std::string& strip(std::string& s, const std::string& chars = blankchars)
        {
            s = s.erase( s.find_last_not_of( chars ) + 1 );
            s = s.erase( 0,s.find_first_not_of( chars ) );
            return s;
        }

        template<typename T>
        struct ex_to_nested_series
        {
            static T convert
            (
                expression_t expression,
                const std::vector<std::string>& expansion_parameters,
                const std::vector<int>& numbers_of_orders,
                const std::vector<int>& offsets = {},
                size_t current_regulator_index = 0,
                bool truncated = false
            )
            {
                assert(current_regulator_index == numbers_of_orders.size());
                assert(current_regulator_index == expansion_parameters.size());
                return expression.begin()->second;
            }
        };

        template<typename T>
        struct ex_to_nested_series<secdecutil::Series<T>>
        {
            static secdecutil::Series<T> convert
            (
                expression_t expression,
                const std::vector<std::string>& expansion_parameters,
                const std::vector<int>& numbers_of_orders,
                const std::vector<int>& offsets = {},
                size_t current_regulator_index = 0,
                bool truncated = false
            )
            {
                int offset = offsets.empty() ? 0 : offsets.at(current_regulator_index);
                std::vector<T> content;
                content.reserve(numbers_of_orders.at(current_regulator_index));
                for (int i = 0; i < numbers_of_orders.at(current_regulator_index); ++i)
                {
                    // Get terms of order i in regulator current_regulator_index (dropping current regulator)
                    powerlist_t lower_bound( expression.begin()->first.size(), 0);
                    powerlist_t upper_bound( expression.begin()->first.size(), 0);
                    lower_bound[0] = i;
                    upper_bound[0] = i+1;
                    expression_t::iterator itlow = expression.lower_bound(lower_bound);
                    expression_t::iterator itup = expression.lower_bound(upper_bound);
                    expression_t subexpression;
                    for (auto it=itlow; it!=itup; ++it)
                    {
                        powerlist_t powerlist = powerlist_t(it->first.begin()+1, it->first.end()); // drop current regulator
                        subexpression[powerlist] = it->second;
                    }
                    if(subexpression.empty())
                        subexpression[powerlist_t( expression.begin()->first.size() - 1, 0)] = "0";
                    content.push_back
                    (
                        ex_to_nested_series<T>::convert
                        (
                            subexpression,
                            expansion_parameters, numbers_of_orders, offsets,
                            current_regulator_index + 1
                        )
                    );
                }
                return {
                            offset,
                            numbers_of_orders.at(current_regulator_index)+offset-1,
                            content, truncated, expansion_parameters.at(current_regulator_index)
                        };
            }
        };

        template<template<typename> class nested_series_t, typename real_t, typename complex_t>
        nested_series_t<complex_t> read_coefficient
        (
            std::istream& stream,
            const std::vector<int>& required_orders,
            const std::vector<std::string>& names_of_regulators,
            const std::vector<std::string>& names_of_real_parameters,
            const std::vector<std::string>& names_of_complex_parameters,
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters
        )
        {
            std::string line_in, line, exprname, exprval, remainder;
            
            Exparse parser;
            
            parser.symbol_table = names_of_regulators;

            parser.substitution_table["I"] = rational_t("0","1");// imaginary unit
            for(int i = 0; i<names_of_real_parameters.size(); i++)
                parser.substitution_table[names_of_real_parameters.at(i)] = rational_t(real_parameters.at(i));
            for(int i = 0; i<names_of_complex_parameters.size(); i++)
                parser.substitution_table[names_of_complex_parameters.at(i)] = rational_t(complex_parameters.at(i).real(),complex_parameters.at(i).imag());
            
            nested_series_t<rational_t> numerator = ex_to_nested_series<nested_series_t<rational_t>>::convert(parser.parse_expression("1"), parser.symbol_table, std::vector<int>(required_orders.size(),1),{},0); // = 1
            nested_series_t<rational_t> denominator = ex_to_nested_series<nested_series_t<rational_t>>::convert(parser.parse_expression("1"), parser.symbol_table, std::vector<int>(required_orders.size(),1),{},0); // = 1

            std::vector<int> offsets;
            offsets.reserve(names_of_regulators.size());
            bool regulator_factor_processed = false;
            bool last_processed = false;
            while (std::getline(stream, line_in, ';'))
            {
                assert( !last_processed );

                secdecutil::exparse::strip(line_in);
                if (line_in.empty()) {
                    last_processed = true;
                    continue; // std::getline should return false which breaks the loop
                }
                
                // replace newline by spaces
                line.clear();
                for (auto& c : line_in)
                {
                    if (c != '\\' && secdecutil::exparse::blankchars.find(c) == std::string::npos)
                        line += c;
                }

                // split by equal sign (=)
                std::istringstream linestream(line);
                std::getline(linestream, exprname, '=');
                secdecutil::exparse::strip(exprname);
                assert(std::getline(linestream, exprval, '='));
                linestream >> remainder;
                assert(remainder.length() == 0);

                //std::cout << "line_in " << line_in << std::endl;
                //std::cout << "exprname " << exprname << std::endl;
                //std::cout << "exprval " << exprval << std::endl;
                
                if(!regulator_factor_processed)
                {
                    // Step 1: Get regulator_factor
                    if(exprname == "regulator_factor")
                    {
                        expression_t regulator_factor_expression = parser.parse_expression(exprval);
                        assert(regulator_factor_expression.size()==1); // Expect only 1 term
                        offsets = std::vector<int>(regulator_factor_expression.begin()->first.begin(), regulator_factor_expression.begin()->first.end());
                        regulator_factor_processed = true;
                        //std::cout << "#regulator_factor_processed: ";
                        //for(auto elem: offsets) std::cout << elem << " ";
                        //std::cout << std::endl;
                        numerator *= ex_to_nested_series<nested_series_t<rational_t>>::convert(parser.parse_expression("1"), parser.symbol_table, std::vector<int>(required_orders.size(),1),offsets,0); // = regulator_factor
                    }
                } else {
                    // Step 2: Process numerators and denominators
                    std::string token;
                    std::istringstream exprval_stream(exprval);
                    while(std::getline(exprval_stream, token, ',')) {
                        if(exprname == "numerator")
                        {
                            numerator *= ex_to_nested_series<nested_series_t<rational_t>>::convert(parser.parse_expression(token), parser.symbol_table, required_orders);
                            //std::cout << "#numerator: " << numerator << std::endl;
                        } else if (exprname == "denominator")
                        {
                            denominator *= ex_to_nested_series<nested_series_t<rational_t>>::convert(parser.parse_expression(token), parser.symbol_table, required_orders);
                            //std::cout << "#denominator: " << denominator << std::endl;
                        } else
                        {
                            throw unknown_expression_error("Encountered unknown expression name while parsing coefficient: " + exprname);
                        }
                    }
                }
            }
            nested_series_t<rational_t> coefficient = numerator/denominator;

            //std::cout << numerator << std::endl;
            //std::cout << denominator << std::endl;
            //std::cout << coefficient << std::endl;
            
            const std::function<complex_t(const rational_t&)> to_complex =
                        [](const rational_t& arg)
                        {
                            return complex_t{mpq_get_d(arg.re), mpq_get_d(arg.im)};
                        };
            return secdecutil::deep_apply(coefficient,to_complex);

        }
    };
};

