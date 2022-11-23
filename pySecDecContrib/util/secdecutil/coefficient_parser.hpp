#include <stdlib.h> // mkstemp
#include <unistd.h> // write

#include <cassert> // assert
#include <cmath> // std::frexp
#include <string> // std::string
#include <vector> // std::vector
#include <limits> // std::numeric_limits
#include <fstream> // std::ifstream
#include <algorithm> // std::erase
#include <array> // std::array
#include <regex>  // std::regex, std::regex_replace

#include <gmp.h> // mpq_get_str

#include <secdecutil/series.hpp> // secdecutil::Series
#include <secdecutil/deep_apply.hpp> // secdecutil::deep_apply

#include <exparse.hpp>

namespace secdecutil
{
    namespace exparse
    {
    
        struct unknown_expression_error : public std::runtime_error { using std::runtime_error::runtime_error; };
    
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
                size_t current_regulator_index = 0
            )
            {
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
                size_t current_regulator_index = 0
            )
            {
                // Print exponent vectors : coefficient
                // for(auto it = expression.cbegin(); it != expression.cend(); ++it)
                // {
                //     for(auto elem : it->first)
                //         std::cout << elem << " ";
                //     std::cout << ": " << it->second << std::endl;
                // }
                int min_order = expression.begin()->first[0];
                int max_order = expression.rbegin()->first[0];
                std::vector<T> content;
                content.reserve(max_order - min_order + 1);
                for (int i = min_order; i <= max_order; ++i)
                {
                    // Get terms of order i in regulator current_regulator_index (dropping current regulator)
                    powerlist_t lower_bound( expression.begin()->first.size(), std::numeric_limits<long long int>::min());
                    powerlist_t upper_bound( expression.begin()->first.size(), std::numeric_limits<long long int>::min());
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
                            expansion_parameters,
                            current_regulator_index + 1
                        )
                    );
                }
                return {min_order, max_order, content, false, expansion_parameters.at(current_regulator_index)};
            }
        };
    
        template<typename real_t, typename complex_t>
        std::string run_ginsh
        (
            const std::string& filename,
            const std::vector<int>& required_orders,
            const std::vector<std::string>& names_of_regulators,
            const std::vector<std::string>& names_of_real_parameters,
            const std::vector<std::string>& names_of_complex_parameters,
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters
        )
        {
            // Get SECDEC_CONTRIB
            std::string secdec_contrib;
            char *secdec_contrib_char = std::getenv("SECDEC_CONTRIB");
            if(secdec_contrib_char == nullptr) {
                secdec_contrib = std::string(SECDEC_CONTRIB);
            } else {
                secdec_contrib = std::string(secdec_contrib_char);
            }

            // Get tmp directory
            std::string tmp_file;
            char *tmp_file_char = std::getenv("TMP");
            if(tmp_file_char == nullptr) {
                tmp_file = std::string("/tmp");
            } else {
                tmp_file = std::string(tmp_file_char);
            }
            tmp_file += "/ginsh_tmp.XXXXXXXXXXXXXXXXXXXX";

            // Create tmp ginsh file
            char *tmp_ginsh_filename = &tmp_file[0]; // template for our file.        
            int fd = mkstemp(&tmp_file[0]);    // Creates and opens a new temp file r/w.
            if(fd<1)
                throw std::runtime_error("failed to open temporary file" + tmp_file + " in coefficient_parser");

            // Populate ginsh file
            std::stringstream ginsh_ss;
            for(const auto& regulator_name : names_of_regulators)
            {
                ginsh_ss << "real_symbols('" << regulator_name << "'):\n";
            }
            for(int i = 0; i < names_of_real_parameters.size(); i++)
                ginsh_ss << names_of_real_parameters.at(i) << "=" << mpq_get_str(NULL,10,rational_t(real_parameters.at(i)).re) << ":\n";
            for(int i = 0; i < names_of_complex_parameters.size(); i++)
                ginsh_ss << names_of_complex_parameters.at(i) << "=" << mpq_get_str(NULL,10,rational_t(complex_parameters.at(i).real(),0).re) << "+I*(" <<  mpq_get_str(NULL,10,rational_t(0,complex_parameters.at(i).imag()).im) << "):\n";
            std::ifstream coefficient_file(filename);
            if (coefficient_file.is_open())
            {
                std::regex e_power("(\\*\\*)");
                std::string line;
                while (std::getline(coefficient_file, line)) {
                    line = std::regex_replace(line, e_power, "^"); // replace ** -> ^ (maintains backwards compatibility)
                    ginsh_ss << "EXPR=" << strip(line) << ":\n";
                }
                coefficient_file.close();
            }
            ginsh_ss << "EXPR2=expand(";
            for(int i = 0; i < names_of_regulators.size(); i++)
            {
                ginsh_ss << "series_to_poly(series(";
            }
            ginsh_ss << "EXPR";
            for(int i = 0; i < names_of_regulators.size(); i++)
            {
                ginsh_ss << "," << names_of_regulators.at(i) << "," << required_orders.at(i)+1 << "))";
            }
            ginsh_ss << "):\n";
            ginsh_ss << "START;" << std::endl;
            ginsh_ss << "expand(real_part(EXPR2)+I_*imag_part(EXPR2));\n";
            ginsh_ss << "quit;\n";
            std::string ginsh_str = ginsh_ss.str();
            ssize_t w = write(fd, ginsh_str.c_str(), ginsh_str.size());
            if (w != ginsh_str.size()) {
                throw std::runtime_error("failed to write to ginsh in coefficient_parser");
            }
            //std::cout << ginsh_str;

            // Run ginsh
            std::string cmd = secdec_contrib + "/bin/ginsh " + std::string(tmp_ginsh_filename);
            std::string result;
            std::array<char, 256> buffer;
            std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
            if (!pipe) {
                // Delete ginsh file
                close(fd);
                unlink(tmp_ginsh_filename);
                throw std::runtime_error("popen() failed when launching ginsh");
            }
            while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
                result += buffer.data();
            }

            // Delete ginsh file
            close(fd);
            unlink(tmp_ginsh_filename);

            // Strip nonsense from ginsh output
            size_t start = result.find("START");
            if (start==std::string::npos)
                throw std::runtime_error("failed to parse ginsh output in coefficient_parser");
            result.erase(0,start+6); // welcome message
            result.pop_back(); // final newline
            result.erase(std::remove(result.begin(), result.end(), '('), result.end()); // (
            result.erase(std::remove(result.begin(), result.end(), ')'), result.end()); // )

            //std::cout << result << std::endl;
            return result;
        }

        template<template<typename> class nested_series_t, typename real_t, typename complex_t>
        nested_series_t<complex_t> read_coefficient
        (
            const std::string& filename,
            const std::vector<int>& required_orders,
            const std::vector<std::string>& names_of_regulators,
            const std::vector<std::string>& names_of_real_parameters,
            const std::vector<std::string>& names_of_complex_parameters,
            const std::vector<real_t>& real_parameters,
            const std::vector<complex_t>& complex_parameters
        )
        {
            Exparse parser;
            parser.symbol_table = names_of_regulators;
            parser.substitution_table["I_"] = rational_t("0","1");// imaginary unit
            for(int i = 0; i<names_of_real_parameters.size(); i++)
                parser.substitution_table[names_of_real_parameters.at(i)] = rational_t(real_parameters.at(i));
            for(int i = 0; i<names_of_complex_parameters.size(); i++)
                parser.substitution_table[names_of_complex_parameters.at(i)] = rational_t(complex_parameters.at(i).real(),complex_parameters.at(i).imag());
            nested_series_t<rational_t> coefficient = ex_to_nested_series<nested_series_t<rational_t>>::convert(
                                                        parser.parse_expression(
                                                                                run_ginsh(
                                                                                            filename,
                                                                                            required_orders,
                                                                                            names_of_regulators,
                                                                                            names_of_real_parameters,
                                                                                            names_of_complex_parameters,
                                                                                            real_parameters,
                                                                                            complex_parameters
                                                                                        )
                                                                            ), 
                                                        names_of_regulators
                                                    );
            
            const std::function<complex_t(const rational_t&)> to_complex =
                        [](const rational_t& arg)
                        {
                            return complex_t{mpq_get_d(arg.re), mpq_get_d(arg.im)};
                        };
            return secdecutil::deep_apply(coefficient,to_complex);
        }
    };
};

