#ifndef EXPARSE_H
#define EXPARSE_H

#include <string> // string
#include <unordered_map> // unordered_map
#include <cassert> // assert
#include <memory> // shared_ptr, make_shared
#include <cstddef> // size_t
#include <cstring> // strcmp
#include <vector> // vector
#include <map> // map
//#include <cstdlib> // atof

#include <gmp.h>

#include <iostream>

// complex<mpq_class>
struct mpqc_class
{
    mpq_t re;
    mpq_t im;
    
    // Comparator Operators
    friend bool operator==(const mpqc_class& lhs, const mpqc_class& rhs){ return (mpq_equal(lhs.re,rhs.re) && mpq_equal(lhs.im,rhs.im)); }
    friend bool operator!=(const mpqc_class& lhs, const mpqc_class& rhs){ return !(lhs == rhs); }
    
    // Unary Operators
    mpqc_class operator-() const {
        mpq_t neg_re;
        mpq_t neg_im;
        mpq_init(neg_re);
        mpq_init(neg_im);
        mpq_neg(neg_re,this->re);
        mpq_neg(neg_im,this->im);
        mpqc_class neg_this(mpq_get_str(NULL,10,neg_re),mpq_get_str(NULL,10,neg_im));
        mpq_clear(neg_re);
        mpq_clear(neg_im);
        return neg_this;
    }
    mpqc_class operator+() const { return *this; }

    // Compound assignment operators
    mpqc_class& operator+=(const mpqc_class& rhs)
    {
        mpq_add(this->re,this->re,rhs.re);
        mpq_add(this->im,this->im,rhs.im);
        return *this;
    }
    mpqc_class& operator-=(const mpqc_class& rhs)
    {
        mpq_sub(this->re,this->re,rhs.re);
        mpq_sub(this->im,this->im,rhs.im);
        return *this;
    }
    mpqc_class& operator*=(const mpqc_class& rhs)
    {
        mpqc_class lhs = *this;
        
        mpq_t nr1;
        mpq_t nr2;
        mpq_t ni1;
        mpq_t ni2;
        mpq_init(nr1);
        mpq_init(nr2);
        mpq_init(ni1);
        mpq_init(ni2);
        
        // re = lhs.re*rhs.re - lhs.im*rhs.im
        mpq_set(nr1,lhs.re);
        mpq_mul(nr1,nr1,rhs.re);
        mpq_set(nr2,lhs.im);
        mpq_mul(nr2,nr2,rhs.im);
        mpq_sub(nr1,nr1,nr2);
        
        // im = lhs.re*rhs.im + lhs.im*rhs.re
        mpq_set(ni1,lhs.re);
        mpq_mul(ni1,ni1,rhs.im);
        mpq_set(ni2,lhs.im);
        mpq_mul(ni2,ni2,rhs.re);
        mpq_add(ni1,ni1,ni2);

        mpq_set(this->re,nr1);
        mpq_set(this->im,ni1);
        
        mpq_clear(nr1);
        mpq_clear(nr2);
        mpq_clear(ni1);
        mpq_clear(ni2);

        return *this;
    }
    mpqc_class& operator/=(const mpqc_class& rhs)
    {
        mpqc_class lhs = *this;
        
        mpq_t nr1;
        mpq_t nr2;
        mpq_t ni1;
        mpq_t ni2;
        mpq_t d1;
        mpq_t d2;
        mpq_init(nr1);
        mpq_init(nr2);
        mpq_init(ni1);
        mpq_init(ni2);
        mpq_init(d1);
        mpq_init(d2);
        
        // re = lhs.re*rhs.re + lhs.im*rhs.im
        mpq_set(nr1,lhs.re);
        mpq_mul(nr1,nr1,rhs.re);
        mpq_set(nr2,lhs.im);
        mpq_mul(nr2,nr2,rhs.im);
        mpq_add(nr1,nr1,nr2);
        
        // im = lhs.re*rhs.im + lhs.im*rhs.re
        mpq_set(ni1,lhs.im);
        mpq_mul(ni1,ni1,rhs.re);
        mpq_set(ni2,lhs.re);
        mpq_mul(ni2,ni2,rhs.im);
        mpq_sub(ni1,ni1,ni2);
        
        // den = rhs.re*rhs.re + rhs.im*rhs.im
        mpq_set(d1,rhs.re);
        mpq_mul(d1,d1,rhs.re);
        mpq_set(d2,rhs.im);
        mpq_mul(d2,d2,rhs.im);
        mpq_add(d1,d1,d2);
        
        mpq_div(nr1,nr1,d1);
        mpq_div(ni1,ni1,d1);
        
        mpq_set(this->re,nr1);
        mpq_set(this->im,ni1);
        
        mpq_clear(nr1);
        mpq_clear(nr2);
        mpq_clear(ni1);
        mpq_clear(ni2);
        mpq_clear(d1);
        mpq_clear(d2);

        return *this;
    }

    // Binary operators
    friend mpqc_class operator+(mpqc_class lhs, const mpqc_class& rhs)
    {
        lhs += rhs;
        return lhs;
    }
    friend mpqc_class operator-(mpqc_class lhs, const mpqc_class& rhs)
    {
        lhs -= rhs;
        return lhs;
    }
    friend mpqc_class operator*(mpqc_class lhs, const mpqc_class& rhs)
    {
        lhs *= rhs;
        return lhs;
    }
    friend mpqc_class operator/(mpqc_class lhs, const mpqc_class& rhs)
    {
        lhs /= rhs;
        return lhs;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const mpqc_class& num)
    {
        os << "(" << mpq_get_str(NULL,10,num.re) << "," << mpq_get_str(NULL,10,num.im) << ")"; // std::complex<> style output
        return os;
    }
    
    mpqc_class& operator= (const mpqc_class& rhs)
    {
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set(this->re,rhs.re);
        mpq_set(this->im,rhs.im);
        return *this;
    };
    
    mpqc_class& operator=(mpqc_class&& rhs)
    {
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set(this->re,rhs.re);
        mpq_set(this->im,rhs.im);
        return *this;
    }

    
    // Constructors
    mpqc_class() {
        mpq_init(this->re);
        mpq_init(this->im);
    };
    mpqc_class(const char* re) {
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set_str(this->re,re,10);
        mpq_canonicalize(this->re);
    };
    mpqc_class(const char* re, const char* im){
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set_str(this->re,re,10);
        mpq_set_str(this->im,im,10);
        mpq_canonicalize(this->re);
        mpq_canonicalize(this->im);
    };
    mpqc_class(double re) {
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set_d(this->re,re);
        mpq_canonicalize(this->re);
    };
    mpqc_class(double re, double im){
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set_d(this->re,re);
        mpq_set_d(this->im,im);
        mpq_canonicalize(this->re);
        mpq_canonicalize(this->im);
    };
    // Copy Constructor
    mpqc_class(const mpqc_class& other)
    {
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set(this->re,other.re);
        mpq_set(this->im,other.im);
    };
    // Move Constructor
    mpqc_class(mpqc_class&& other)
    {
        mpq_init(this->re);
        mpq_init(this->im);
        mpq_set(this->re,other.re);
        mpq_set(this->im,other.im);
    }
    // Destructor
    ~mpqc_class()
    {
        mpq_clear(this->re);
        mpq_clear(this->im);
    };
};

class Exparse
{
private:
    
    typedef mpqc_class rational_t;
    typedef long long int int_t;
    
    enum Operation { add, subtract, multiply, divide, inverse};
    
    const rational_t rational_minus_one  = mpqc_class("-1");
    const rational_t rational_zero = mpqc_class("0");
    const rational_t rational_one = mpqc_class("1");
    const rational_t rational_two = mpqc_class("2");
    const int_t integer_zero = 0;
    const int_t integer_one = 1;
    const int_t integer_two = 2;
    
    rational_t term_buffer = rational_zero;
    rational_t symbol_buffer = rational_zero;
    rational_t number_buffer = rational_zero;
    rational_t pow_buffer = rational_zero;
    
    int_t int_buffer = integer_one;

    std::vector<int_t> symbol_orders_buffer;
    std::size_t symbol_order_buffer;
    bool symbol_order_altered = false;
    
    struct slice_t
    {
        std::shared_ptr<std::string> expression;
        std::size_t pos;
        std::size_t len;
        
        const char& operator[](const std::size_t index) const { return (*expression)[pos+index]; }
        
        friend inline bool operator==(const slice_t& lhs, const std::string& rhs)
        {
            if (lhs.len != rhs.length())
            return false;
            for (std::size_t i=0; i<lhs.len; ++i)
            {
                if( lhs[i] != rhs[i] )
                return false;
            }
            return true;
        }
        friend inline bool operator!=(const slice_t& lhs, const std::string& rhs){ return !(lhs == rhs); }
        friend std::ostream& operator<<(std::ostream& os, const slice_t& slice) { os << slice.expression->substr(slice.pos,slice.len); return os; }
    };
    
    void exparse_apply(rational_t& result, const rational_t& symbol, const Operation op)
    {
        if(symbol_order_altered)
        {
            symbol_order_altered = false;
            symbol_orders_buffer[symbol_order_buffer] += 1; // exparse_apply is called when no power is set
            return;
        }
        switch(op)
        {
            case add: result *= symbol; break;
            case subtract: result *= -symbol; break;
            case multiply: result *= symbol; break;
            case divide: result /= symbol; break;
            case inverse: result = symbol/result; break;
        }
    }
    
    void exparse_pow(rational_t& base, const int_t& exponent, rational_t& result, Operation op)
    {
        if(symbol_order_altered)
        {
            symbol_order_altered = false;
            symbol_orders_buffer[symbol_order_buffer] += exponent;
            return;
        }

        if (exponent < integer_zero)
        {
            exparse_apply(result, rational_one, inverse); // 1/result
            exparse_pow(base, -exponent, result, op);
            exparse_apply(result, rational_one, inverse); // 1/result
            return;
        }
        else if (exponent == integer_zero)
        {
            return; // nothing to do
        }
        else if (exponent == integer_one)
        {
            exparse_apply(result, base, op);
            return;
        }
        else if (exponent == integer_two)
        {
            exparse_apply(result, base*base, op);
            return;
        }
        
        // Set power buffer to one (all symbols will be multiplied on to power buffer)
        pow_buffer = rational_one;
        
        int_t tmp_exponent = exponent;
        while(tmp_exponent > 0)
        {
            if(tmp_exponent & 1) // exponent is odd
            pow_buffer *= base;
            base *= base;
            tmp_exponent = tmp_exponent >> 1;
        }
        exparse_apply(result, pow_buffer, op);
    }
    
    void to_number(slice_t& symbol, rational_t& result)
    {
        char store;
        if ( symbol.pos+symbol.len != symbol.expression->length() )
        {
            // Store character after symbol and replace with '\0'
            store =(*symbol.expression)[symbol.pos+symbol.len];
            (*symbol.expression)[symbol.pos+symbol.len] = '\0';
        }
        
        // Note: would be nicer to use find() but does not work correctly for pointer types
        bool found = false;
        // Check if symbol is in substitution table
        for( const std::pair<std::string,rational_t> substitution_table_element: substitution_table)
        {
            if ( strcmp(substitution_table_element.first.c_str(), symbol.expression->c_str()+symbol.pos) == 0)
            {
                result = substitution_table_element.second;
                found = true;
                break;
            }
        }
        // Check if symbol is in symbol_table
        if(!found)
        {
            for( std::size_t i = 0; i < symbol_table.size(); i++)
            {
                if ( strcmp(symbol_table[i].c_str(), symbol.expression->c_str()+symbol.pos) == 0)
                {
                    symbol_order_altered = true;
                    symbol_order_buffer = i;
                    found = true;
                    break;
                }
            }
        }

        if(!found)
        {
            // Parse symbol as rational
            result = rational_t(symbol.expression->c_str()+symbol.pos);
        }

        if ( symbol.pos+symbol.len != symbol.expression->length() )
        {
            // Restore character after symbol
            (*symbol.expression)[symbol.pos+symbol.len] = store;
        }
    }
    
    void to_int(slice_t& symbol, int_t& result)
    {
        if ( symbol.pos+symbol.len != symbol.expression->length() )
        {
            const char store =(*symbol.expression)[symbol.pos+symbol.len];
            (*symbol.expression)[symbol.pos+symbol.len] = '\0';
            result = std::atoll(symbol.expression->c_str()+symbol.pos);
            (*symbol.expression)[symbol.pos+symbol.len] = store;
        }
        else
        {
            result = std::atoll(symbol.expression->c_str()+symbol.pos);
        }
    }
    
    void parse_symbol(slice_t& symbol, rational_t& result)
    {
        // std::size_t symbol_pos = symbol.pos;
        // std::size_t symbol_len = symbol.len;

        // Parse operator
        Operation op = add;
        if( symbol[0] == '+')
        {
            symbol.pos++;
            symbol.len--;
            op = add;
        }
        else if( symbol[0] == '-')
        {
            symbol.pos++;
            symbol.len--;
            op = subtract;
        }
        else if( symbol[0] == '*')
        {
            symbol.pos++;
            symbol.len--;
            op = multiply;
        }
        else if( symbol[0] == '/')
        {
            symbol.pos++;
            symbol.len--;
            op = divide;
        }
        
        // Parse power
        for (std::size_t i=0; i<symbol.len; ++i)
        {
            if( symbol[i] == '^')
            {
                slice_t base_slice = {symbol.expression,symbol.pos,i};
                slice_t exponent_slice = {symbol.expression,symbol.pos+i+1,symbol.len-(i+1)};
                
                to_number(base_slice,number_buffer);
                to_int(exponent_slice,int_buffer);
                exparse_pow(number_buffer, int_buffer, result, op);
                
                return;
            }
        }
        
        // Parse non-power
        to_number(symbol,number_buffer);
        exparse_apply(result,number_buffer,op);
        return;
    }
    
    void parse_term(slice_t& term)
    {
        std::size_t term_pos = term.pos;
        std::size_t term_len = term.len;
        
        // Set term to one (all symbols will be multiplied on to term)
        term_buffer = rational_one;

        // Reset symbol_orders_buffer
        symbol_orders_buffer.assign(symbol_table.size() ,0);
        
        // Parse term
        std::size_t reading_point = 0;
        for (std::size_t i=0; i<term.len; ++i)
        {
            if( term[i] == '*' || term[i] == '/')
            {
                // Parse term
                term.pos = term_pos+ reading_point;
                term.len = i-reading_point;
                parse_symbol(term, term_buffer);
                
                // Update reading point
                reading_point = i;
                
                // Reset term to original parameters
                term.pos = term_pos;
                term.len = term_len;
            }
            else if ( i == term.len-1 )
            {
                // Parse term
                term.pos = term_pos+ reading_point;
                term.len = i+1-reading_point;
                parse_symbol(term, term_buffer);
                
                // Update reading point
                reading_point = i+1;
                
                // Reset term to original parameters
                term.pos = term_pos;
                term.len = term_len;
            }
        }
    }

    void add_term(std::map<std::vector<int_t>,rational_t>& result)
    {
        std::map<std::vector<int_t>,rational_t>::iterator lb = result.lower_bound(symbol_orders_buffer);
        if(lb != result.end() && !result.key_comp()(symbol_orders_buffer, lb->first))
        {
            // symbol_orders exists in result, add term to result
            lb->second += term_buffer;
        }
        else
        {
            // symbol_orders does not exist in result, set term equal to term
            result.insert(lb, std::map<std::vector<int_t>,rational_t>::value_type(symbol_orders_buffer, term_buffer));
        }

        // Reset buffers
    }

    void parse_line(slice_t& line, std::map<std::vector<int_t>,rational_t>& result)
    {
        std::size_t line_pos = line.pos;
        std::size_t line_len = line.len;

        symbol_orders_buffer.assign(symbol_table.size(), 0); // initialise symbol_orders_buffer

        // Parse line
        std::size_t reading_point = 0;
        for (std::size_t i=0; i<line.len; ++i)
        {
            if( i != 0 && (line[i] == '+' || (line[i] == '-' && line[i-1] != '^') ))
            {
                // Parse term
                line.pos = line_pos + reading_point;
                line.len = i-reading_point;

                parse_term(line);
                add_term(result);

                // Update reading point
                reading_point = i;
                
                // Reset line to original parameters
                line.pos = line_pos;
                line.len = line_len;
            }
            else if (i == line.len-1)
            {
                // Parse term
                line.pos = line_pos + reading_point;
                line.len = i+1-reading_point;
                
                parse_term(line);
                add_term(result);

                // Update reading point
                reading_point = i+1;
                
                // Reset line to original parameters
                line.pos = line_pos;
                line.len = line_len;
            }
        }
    }
    
    void parse_sanity(const slice_t& line)
    {
        for (std::size_t i=0; i<line.expression->length(); ++i)
        {
            assert(line[i]!='(');
            assert(line[i]!=')');
            assert(line[i]!=' ');
            assert(line[i]!=';');
        }
        assert(line[line.expression->length()-1] != '+');
        assert(line[line.expression->length()-1] != '-');
        assert(line[line.expression->length()-1] != '*');
        assert(line[line.expression->length()-1] != '^');
    }
    
public:

    std::vector<std::string> symbol_table;
    std::unordered_map<std::string,rational_t> substitution_table;

    std::map<std::vector<int_t>, rational_t> parse_expression(const std::string& expression)
    {
        std::map<std::vector<int_t>, rational_t> result;
        
        slice_t line;
        line.expression = std::make_shared<std::string>(expression);
        line.pos = 0;
        line.len = line.expression->length();
        
        parse_line(line, result);
        
        return result;
    }
    
    // Constructor
    Exparse()
    {
        // Allow buffer sizes to be set?
    }
    
};

#endif

