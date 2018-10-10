#include <cassert> // assert
#include <cmath> // std::frexp
#include <istream> // std::istream
#include <sstream> // std::stringstream
#include <string> // std::string
#include <vector> // std::vector

#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include <ginac/ginac.h>

namespace secdecutil
{
    namespace ginac
    {
        namespace g = GiNaC;

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
                g::ex expression,
                const std::vector<g::symbol>& expansion_parameters,
                const std::vector<int>& numbers_of_orders,
                const std::vector<int>& offsets = {},
                size_t current_regulator_index = 0
            )
            {
                assert(current_regulator_index == numbers_of_orders.size());
                assert(current_regulator_index == expansion_parameters.size());
                return expression;
            }
        };
        template<typename T>
        struct ex_to_nested_series<secdecutil::Series<T>>
        {
            static secdecutil::Series<T> convert
            (
                g::ex expression,
                const std::vector<g::symbol>& expansion_parameters,
                const std::vector<int>& numbers_of_orders,
                const std::vector<int>& offsets = {},
                size_t current_regulator_index = 0
            )
            {
                int offset = offsets.empty() ? 0 : offsets.at(current_regulator_index);
                std::vector<T> content;
                content.reserve(numbers_of_orders.at(current_regulator_index));
                for (int i = 0; i < numbers_of_orders.at(current_regulator_index); ++i)
                {
                    content.push_back
                    (
                        ex_to_nested_series<T>::convert
                        (
                            expression.coeff(expansion_parameters.at(current_regulator_index),i),
                            expansion_parameters, numbers_of_orders, offsets,
                            current_regulator_index + 1
                        )
                    );
                }
                return {
                            offset,
                            numbers_of_orders.at(current_regulator_index)+offset-1,
                            content, true, expansion_parameters.at(current_regulator_index).get_name()
                        };
            }
        };

        inline g::numeric pow_to_unsigned(const g::numeric& base, unsigned exponent)
        {
            if(exponent == 0)
                return 1;

            unsigned half_exponent = exponent / 2;
            g::numeric out = pow_to_unsigned(base, half_exponent);

            out *= out;
            if (2 * half_exponent == exponent) // exponent is even
                return out;
            else // exponent is odd --> need another factor of the base due to integer division above
                return out * base;

        }

        inline g::numeric rationalize(const double input)
        {
            const g::numeric radix(2);
            const long long radix_to_signficand_digits = 9007199254740992; // 2^53
            int exp;
            double significand = frexp(input,&exp);

            g::numeric converted_numerator = g::numeric(std::to_string(static_cast<long long>(significand * radix_to_signficand_digits)).c_str());
            g::numeric converted_denominator = g::numeric(std::to_string(radix_to_signficand_digits).c_str());
            if (exp > 0)
                converted_numerator *= pow_to_unsigned(radix,exp);
            else if (exp < 0)
                converted_denominator *= pow_to_unsigned(radix,-exp);
            // else exp == 0, nothing to do

            assert(input - (converted_numerator/converted_denominator) <= std::numeric_limits<double>::epsilon()); // check rationalization

            return converted_numerator/converted_denominator;
        }

        template<template<typename> class nested_series_t, typename real_t, typename complex_t>
        nested_series_t<complex_t> read_coefficient // TODO: test cases and documentation
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

            // make a symtab with the regulators, and with the real and complex parameters
            g::symtab expressions;

            // convert the regulator names to ginac symbols
            std::vector<g::symbol> regulator_symbols;

            // populate "expressions" and "regulator_symbols"
            regulator_symbols.reserve(names_of_regulators.size());
            for(const std::string& name : names_of_regulators)
            {
                regulator_symbols.push_back(g::symbol(name.c_str()));
                expressions[name.c_str()] = regulator_symbols.back();
            }
            for(unsigned int i = 0; i < names_of_real_parameters.size(); ++i)
                expressions[names_of_real_parameters.at(i).c_str()] = secdecutil::ginac::rationalize(real_parameters.at(i));
            for(unsigned int i = 0; i < names_of_complex_parameters.size(); ++i)
            {
                complex_t parameter = complex_parameters.at(i);
                expressions[names_of_complex_parameters.at(i).c_str()] =
                    rationalize(parameter.real()) + rationalize(parameter.imag()) * g::I;
            }

            bool last_processed = false;
            while (std::getline(stream, line_in, ';'))
            {
                assert( !last_processed );

                secdecutil::ginac::strip(line_in);
                if (line_in.empty()) {
                    last_processed = true;
                    continue; // std::getline should return false which breaks the loop
                }

                // replace newline by spaces
                line.clear();
                for (auto& c : line_in)
                {
                    if (c != '\\' && secdecutil::ginac::blankchars.find(c) == std::string::npos)
                        line += c;
                }

                // split by equal sign (=)
                std::istringstream linestream(line);
                std::getline(linestream, exprname, '='); secdecutil::ginac::strip(exprname);
                assert(std::getline(linestream, exprval, '='));
                linestream >> remainder; assert(remainder.length() == 0);

                // parse expr
                g::parser reader(expressions);
                expressions[exprname] = reader(exprval).expand();
            }


            // account for "regulator_factor" in "numerator"
            g::ex regulator_factor = expressions["regulator_factor"];
            std::vector<int> offsets; offsets.reserve(names_of_regulators.size());
            for(const g::symbol& regulator : regulator_symbols)
            {
                int degree = regulator_factor.ldegree(regulator);
                assert(degree == regulator_factor.degree(regulator));
                offsets.push_back(degree);
            }

            // convert numerator and denominator to nested_series_t<g::ex>
            nested_series_t<g::ex> numerator =
                secdecutil::ginac::ex_to_nested_series<nested_series_t<g::ex>>::convert
                (
                    expressions["numerator"],
                    regulator_symbols,
                    required_orders,
                    offsets
                );
            nested_series_t<g::ex> denominator =
                secdecutil::ginac::ex_to_nested_series<nested_series_t<g::ex>>::convert
                (
                    expressions["denominator"],
                    regulator_symbols,
                    required_orders
                );

            // convert the g::ex to complex_t
            const std::function<complex_t(const g::ex&)> ex_to_complex =
            [](const g::ex& arg)
            {
                assert( g::is_a<g::numeric>(arg) ); // "ex_to" without assertion is unsafe
                const auto& value = g::ex_to<g::numeric>(arg);
                return complex_t{value.real().to_double(), value.imag().to_double()};
            };
            return secdecutil::deep_apply(numerator/denominator, ex_to_complex);
        }

    };
};

