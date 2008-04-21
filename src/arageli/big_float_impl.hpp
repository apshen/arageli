/*****************************************************************************

    big_float_impl.hpp

    This file is a part of the Arageli library.

    Copyright (C) 2005--2006 Alexander Pshenichnikov
    Copyright (C) 2005--2006 Nikolay Santalov
    University of Nizhni Novgorod, Russia

    The Arageli Library is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License version 2
    as published by the Free Software Foundation.

    The Arageli Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.

    We are also open for dual licensing for the whole library or
    for its particular part. If you are interested to get the library
    in this way, i.e. not under the GNU General Public License,
    please contact Arageli Support Service support.arageli@gmail.com.

*****************************************************************************/

/**
    \file big_float_impl.hpp
    \brief Big float number class implementation.

    This module implements a class big_float for representing
    Big Float Numbers, i.e. unlimited precision real numbers with
    floating point.

    Implementation of a Big Float Number.
    An instance of the data type big_float_impl is
    a number of the form s*2^e where s is
    the significant and e is the exponent.
    Both are instances of type big_int.
    There are the special big_float values
    NaN (not a number), pZero, nZero (+0, -0),
    pInf, mInf (+infinity, -infinity).

    Arithmetic on big_float numbers uses two
    parameters: precision and rounding mode
*/

#include "config.hpp"

#include <cmath>
#include "std_import.hpp"
#include "bigar.hpp"
#include "big_int.hpp"

#ifndef _ARAGELI_big_float_impl_h_
#define _ARAGELI_big_float_impl_h_

#undef NAN

#pragma warning ( disable : 4018 )

namespace Arageli
{

/** \name Rounding Modes for big_float*/
//@{

// WARNING! Names are too simple and general.
// TODO: Rename or incapsulate them!
// WARNING! How does it combine with ISO C++ 18.2.1.3?

const int EXACT = 1; ///< Rounds to the closest value
const int ROUND = 2; ///< Uses rules of `'manual'' computations
const int TO_NEAREST = 3; ///< Rounds towards the nearest number
const int TO_ZERO = 4; ///< Rounds towards zero
const int TO_P_INF = 5; ///< Rounds towards +infinity
const int TO_M_INF = 6; ///< Rounds towards -infinity
const int TO_INF = 7;   ///< Rounds outwards zero
//@}

const int PREC_MIN = 2;
const unsigned long PREC_MAX = _Internal::max_digit;
const unsigned long D_PREC_MAX =
    ( long ) (_Internal::max_digit * log ( 2.0l ) / log ( 10.0l ));

class big_float_impl
{

public:

    /**@name Initializations from built-in types. */
    //@{
    big_float_impl (char i)
    {
        from_native_int(i);
    }

    big_float_impl (signed char i)
    {
        from_native_int(i);
    }

    big_float_impl (unsigned char i)
    {
        from_native_int(i);
    }

    big_float_impl (signed short i)
    {
        from_native_int(i);
    }

    big_float_impl (unsigned short i)
    {
        from_native_int(i);
    }

    big_float_impl (signed int i)
    {
        from_native_int(i);
    }

    big_float_impl (unsigned int i)
    {
        from_native_int(i);
    }

    big_float_impl (signed long int i)
    {
        from_native_int(i);
    }

    big_float_impl (unsigned long int i)
    {
        from_native_int(i);
    }

    big_float_impl (bool i)
    {
        from_native_int(i);
    }

    big_float_impl (const big_int &i);

    #ifdef ARAGELI_INT64_SUPPORT
        big_float_impl (signed __int64 i)
        {
            from_native_int(i);
        }

        big_float_impl (unsigned __int64 i)
        {
            from_native_int(i);
        }
    #endif

    #ifdef ARAGELI_LONG_LONG_SUPPORT
        big_float_impl (signed long long i)
        {
            from_native_int(i);
        }

        big_float_impl (unsigned long long i)
        {
            from_native_int(i);
        }
    #endif

    big_float_impl (float f)
    {
        from_native_float(f);
    }

    big_float_impl (double f)
    {
        from_native_float(f);
    }

    big_float_impl (long double f)
    {
        from_native_float(f);
    }

    //@}

    /**@name Other constructors */
    //@{
    big_float_impl();                                   ///< Default constructor
    big_float_impl ( long prec, long mode );            ///< Sets mode to mode and prec to prec
    big_float_impl ( const char *str );                 ///< Converts str to a big_float_impl
    big_float_impl ( const char* str, long p );
    big_float_impl ( const big_float_impl & b );             ///< Makes a copy of a number

    /// Sets mantissa to s and exponenta to e
    big_float_impl
    (
        const big_int &s,
        const big_int &e,
        long mode
    );
    //@}

    ~big_float_impl();

    /**@name Assignment*/
    //@{
    big_float_impl& operator= (const big_int & b);

    big_float_impl& operator= (const char* str)
    {
        return *this = big_float_impl (str,prec);
    }

    big_float_impl& operator= (char i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (signed char i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (unsigned char i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (signed short i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (unsigned short i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (signed int i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (unsigned int i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (signed long int i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (unsigned long int i)
    {
        from_native_int(i);
        return *this;
    }

    big_float_impl& operator= (bool i)
    {
        from_native_int(i);
        return *this;
    }


    #ifdef ARAGELI_INT64_SUPPORT
        big_float_impl& operator= (signed __int64 i)
        {
            from_native_int(i);
            return *this;
        }

        big_float_impl& operator= (unsigned __int64 i)
        {
            from_native_int(i);
            return *this;
        }
    #endif

    #ifdef ARAGELI_LONG_LONG_SUPPORT
        big_float_impl& operator= (signed long long i)
        {
            from_native_int(i);
            return *this;
        }

        big_float_impl& operator= (unsigned long long i)
        {
            from_native_int(i);
            return *this;
        }
    #endif

    big_float_impl& operator= (float f)
    {
        from_native_float(f); return *this;
    }

    big_float_impl& operator= (double f)
    {
        from_native_float(f);
        return *this;
    }

    big_float_impl& operator= (long double f)
    {
        from_native_float(f);
        return *this;
    }

    //@}

    /**@name Convertions to built-in types. */
    //@{

    operator char () const
    {
        return to_native_int<char>();
    }

    operator signed char () const
    {
        return to_native_int<signed char>();
    }

    operator unsigned char () const
    {
        return to_native_int<unsigned char>();
    }

    operator signed short () const
    {
        return to_native_int<signed short>();
    }

    operator unsigned short () const
    {
        return to_native_int<unsigned short>();
    }

    operator signed int () const
    {
        return to_native_int<int>();
    }

    operator unsigned int () const
    {
        return to_native_int<unsigned int>();
    }

    operator signed long int () const
    {
        return to_native_int<signed long>();
    }

    operator unsigned long int () const
    {
        return to_native_int<unsigned long>();
    }

    operator bool () const
    {
        return !is_zero();
    }

    #ifdef ARAGELI_INT64_SUPPORT

    operator signed __int64 () const
    {
        return to_native_int<signed __int64>();
    }

    operator unsigned __int64 () const
    {
        return to_native_int<unsigned __int64>();
    }

    #endif

    #ifdef ARAGELI_LONG_LONG_SUPPORT

        operator signed long long () const
        {
            return to_native_int<signed long long>();
        }

        operator unsigned long long () const
        {
            return to_native_int<unsigned long long>();
        }

    #endif

        operator float () const
        {
            return to_native_float<float>();
        }

        operator double () const
        {
            return to_native_float<double>();
        }

        operator long double () const
        {
            return to_native_float<long double>();
        }

    //@}

    /// Returns the exponent
    big_int get_exponent() const
    {
        return e;
    }

    /// Returns the significant
    big_int get_significant() const
    {
        return s;
    }

    /// Returns signum
    int sign () const;

    /// Swaps two numbers.
    void swap (big_float_impl& bf)
    {
        s.swap(bf.s);
        e.swap(bf.e);
        std::swap(prec, bf.prec);
        std::swap(mode, bf.mode);
    }

    /// Reads a number
    friend std::ostream & operator << (std::ostream & s, const big_float_impl & x);
    /// Writes a number
    friend std::istream & operator >> (std::istream & s, big_float_impl & x);
    /// Formated output
    void out ( std::ostream & os, char c = 'd' ) const ;
    /// Formated input
    void in ( std::istream &ins, char c = 'd' );

    /// Sets the global precision to p
    static void set_global_precision(long p);
    /// Sets the global rounding mode
    static void set_global_round_mode (int m = TO_NEAREST);

    /// Returns bits precision required to guarantee dec precision indicated
    static std::size_t dec2bits_precision ( std::size_t dec_precision )
    {
        return (std::size_t) ((long double)(dec_precision) * log(10.0l) / log(2.0l) + 1.0);
    }

    /// Returns dec precision guaranteed by indicated bits precision
    static std::size_t bits2dec_precision ( std::size_t bits_precision )
    {
        return (std::size_t) ((long double)(bits_precision) * log(2.0l) / log(10.0l));
    }

    /// Sets precision to p
    void set_precision ( long p );
    /// Sets rounding mode
    void set_round_mode ( long m );

    /// Get precision
    long get_precision ( void )
    {
        return prec;
    }

    /// Get rounding mode
    long get_round_mode ( void )
    {
        return mode;
    }

    /// Compares two numbers
    /**
        Returns
        -  0  if a = b,
        - -1  if a < b,
        - +1  if a > b
    */
    friend int cmp(const big_float_impl & a, const big_float_impl & b);
    /// Test for equality
    friend bool operator ==(const big_float_impl & a, const big_float_impl & b);
    /// Test for inequality
    friend bool operator !=(const big_float_impl & a, const big_float_impl & b);
    /// Test for greater
    friend bool operator > (const big_float_impl & a, const big_float_impl & b);
    /// Test for greater than or equal to
    friend bool operator >=(const big_float_impl & a, const big_float_impl & b);
    /// Test for less
    friend bool operator < (const big_float_impl & a, const big_float_impl & b);
    /// Test for less than or equal to
    friend bool operator <=(const big_float_impl & a, const big_float_impl & b);

    friend big_float_impl operator + (const big_float_impl & a);       ///< Unary plus
    friend big_float_impl operator - (const big_float_impl & a);       ///< Unary minus

    /** Addition. Up to prec binary digits in rounding mode mode.
        The parameters prec and mode are optional and have the global defult
        values which can be set by set_global_precision and set_rounding_mode.*/
    friend big_float_impl add(const big_float_impl & b, const big_float_impl & c, long prec, int mode);
    friend big_float_impl add(const big_float_impl & b, const big_float_impl & c, long prec);
    friend big_float_impl add(const big_float_impl & b, const big_float_impl & c);

    /** Subtraction. Up to prec binary digits in rounding mode mode
        The parameters prec and mode are optional and have the global defult
        values which can be set by set_global_precision and set_rounding_mode.*/
    friend big_float_impl sub(const big_float_impl & b, const big_float_impl & c, long prec, int mode);
    friend big_float_impl sub(const big_float_impl & b, const big_float_impl & c, long prec);
    friend big_float_impl sub(const big_float_impl & b, const big_float_impl & c);

    /** Multiplication. Up to prec binary digits in rounding mode mode
        The parameters prec and mode are optional and have the global defult
        values which can be set by set_global_precision and set_rounding_mode.*/
    friend big_float_impl mul(const big_float_impl & b, const big_float_impl & c, long prec, int mode);
    friend big_float_impl mul(const big_float_impl & b, const big_float_impl & c, long prec);
    friend big_float_impl mul(const big_float_impl & b, const big_float_impl & c);

    /** Divizion. Up to prec binary digits
        The parameters prec and mode are optional and have the global defult
        values which can be set by set_global_precision */
    friend big_float_impl div(const big_float_impl & b, const big_float_impl & c, long prec, int mode);
    friend big_float_impl div(const big_float_impl & b, const big_float_impl & c, long prec);
    friend big_float_impl div(const big_float_impl & b, const big_float_impl & c);
    friend big_float_impl divnu(const big_float_impl & b, const big_float_impl & c, long prec, int mode);
    friend big_float_impl div_i ( const big_int & c );
    /** Square rooting. Up to prec binary digits
        The parameters prec and mode are optional and have the global defult
        values which can be set by set_global_precision.*/
    friend big_float_impl fsqrt(const big_float_impl & b,/* const big_float_impl & c,*/ long prec, int mode);
    friend big_float_impl fsqrt(const big_float_impl & b,/* const big_float_impl & c,*/ long prec);
    friend big_float_impl fsqrt(const big_float_impl & b/*, const big_float_impl & c*/);
    friend big_float_impl nfsqrt(const big_float_impl & b,/* const big_float_impl & c,*/ long prec, int mode);

    /// Binary plus. It is equivalent to add(b, c)
    friend big_float_impl operator + (const big_float_impl & b, const big_float_impl & c);
    /// Binary minus. It is equivalent to sub(b, c)
    friend big_float_impl operator - (const big_float_impl & b, const big_float_impl & c);
    /// Multiplication. It is equivalent to mul(b, c)
    friend big_float_impl operator * (const big_float_impl & b, const big_float_impl & c);
    /// Divizion. It is equivalent to div(b, c)
    friend big_float_impl operator / (const big_float_impl & b, const big_float_impl & c);
    /// Combined assignment-addition operator
    friend big_float_impl & operator += (big_float_impl & b, const big_float_impl & c);
    /// Combined assignment-subtraction operator
    friend big_float_impl & operator -= (big_float_impl & b, const big_float_impl & c);
    /// Combined assignment-multiplication operator
    friend big_float_impl & operator *= (big_float_impl & b, const big_float_impl & c);
    /// Combined assignment-division operator
    friend big_float_impl & operator /= (big_float_impl & b, const big_float_impl & c);

    /// Returns the next bigger integer
    friend big_float_impl ceil  (const big_float_impl & a);
    /// Returns the next smaller integer
    friend big_float_impl floor (const big_float_impl & a);
    /// Returns the fractal part of the number
    friend big_float_impl frac  (const big_float_impl & a);
    /// Returns the next bigger integer
    friend big_int iceil (const big_float_impl & a);
    /// Returns the next smaller integer
    friend big_int ifloor(const big_float_impl & a);

    /// Returns random number in the interval 0 <=x < 1 with prec lenght of mantissa
    friend big_float_impl frandom  (long bits);
    /// Returns random number in the interval 0 <=x < 1 with prec lenght of mantissa
    friend big_float_impl frandom  ( );
    /// Returns random number with prec lenght of mantissa and with exp as exponent
    friend big_float_impl frandom1 (long bits, const big_int & exp = 0 );

    /// Returns 1 iff the number is Nan, Infinity or 0
    bool is_special ( void ) const
    {
        return special != FINITE;
    }

    /// Returns 1 iff the number is NaN
    bool is_nan     ( void ) const
    {
        return special == NAN;
    }

    /// Returns 1 iff the number is +Infinity
    bool is_pinf    ( void ) const
    {
        return special == P_INF;
    }

    /// Returns 1 iff the number is -Infinity
    bool is_minf    ( void ) const
    {
        return special == M_INF;
    }

    /// Returns 1 iff the number is +Infinity or -Infinity
    bool is_inf     ( void ) const
    {
        return special == P_INF || special == M_INF;
    }

    /// Returns 1 iff the number is +0
    bool is_pzero   ( void ) const
    {
        return special == P_ZERO;
    }

    /// Returns 1 iff the number is -0
    bool is_mzero   ( void ) const
    {
        return special == M_ZERO;
    }

    /// Returns 1 iff the number is +0
    bool is_zero    ( void ) const
    {
        return special == P_ZERO || special == M_ZERO;
    }

    #if 0
    /** \name Special big_float_impl numbers */
    //@{

        static big_float_impl nan ( void )
        {
            return big_float_impl ("nan");
        }

        static big_float_impl inf ( void )
        {
            return big_float_impl ("inf");
        }

        static big_float_impl m_inf ( void )
        {
            return big_float_impl ("-inf");
        }

        static big_float_impl zero ( void )
        {
            return big_float_impl ( );
        }

        static big_float_impl m_zero ( void )
        {
            return big_float_impl ("-0");
        }

    //@}
    #endif

    typedef enum special_numbers_en
    {
        FINITE = 0,
        NAN,
        P_INF,
        M_INF,
        P_ZERO,
        M_ZERO
    } special_numbers;

    /*
    friend special_numbers operator + ( const special_numbers &sa, const special_numbers &sb );
    friend special_numbers operator - ( const special_numbers &sa, const special_numbers &sb );
    friend special_numbers operator * ( const special_numbers &sa, const special_numbers &sb );
    friend special_numbers operator / ( const special_numbers &sa, const special_numbers &sb );

    friend std::ostream & operator << ( std::ostream &os, const special_numbers &sn );
    friend std::istream & operator >> ( std::istream &is, const special_numbers &sn );
    */

private:

    void normalize_1 ( long prec = big_float_impl::global_prec, long mode = big_float_impl::global_mode );

    // the following type is used inside the big_float_impl unit and implements
    // the storage for a Big Float Number

    special_numbers special;

    big_int e; ///< exponent
    big_int s; //< significant

    long prec; ///< bits precision for this big_float_impl number
    long mode; ///< big_float_impl rounding mode for this number

    static long global_prec; ///< Global bits precision for big_float_impl numbers
    static long global_mode; ///< Global big_float_impl rounding mode

private:

    //  static const std::string get_lower_string ( const special_numbers sn );
    //  static const std::string get_upper_string ( const special_numbers sn );

    template <typename T>
    void from_native_float (const T &f);

    template <typename T>
    void from_native_int (const T &i);

    template <typename T>
    T to_native_int () const;

    template <typename T>
    T to_native_float () const;
};


inline big_float_impl abs (const big_float_impl& x)
{
    return is_negative(x) ? -x : x;
}


template <typename Outiter>
inline void generate_range_helper (big_float_impl& t1, const big_float_impl& t2, Outiter outiter)
{
    generate_range_helper_wo_inc(t1, t2, outiter);
}


inline big_float_impl & operator+= (big_float_impl & b, const big_float_impl & c)
{
    return b = b + c;
}

inline big_float_impl & operator-= (big_float_impl & b, const big_float_impl & c)
{
    return b = b - c;
}

inline big_float_impl & operator*= (big_float_impl & b, const big_float_impl & c)
{
    return b = b * c;
}

inline big_float_impl & operator/= (big_float_impl & b, const big_float_impl & c)
{
    return b = b / c;
}


big_float_impl fsqrt( const big_float_impl & b );


}

namespace std
{

inline Arageli::big_float_impl abs(const Arageli::big_float_impl &x)
{
    return Arageli::abs(x);
}

inline Arageli::big_float_impl sqrt(const Arageli::big_float_impl &x)
{
    return Arageli::fsqrt(x);
}

}

#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_BIG_FLOAT
    #include "big_float_impl.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_BIG_FLOAT
#endif

#endif