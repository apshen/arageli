/*****************************************************************************

    big_float_impl.cpp

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
    \file big_float_impl.cpp
    \brief The big_float_impl.hpp file stuff implementation.

    Big Float Numbers implementation.
*/


#include "config.hpp"

#define CHECK_PREC(p)  {}  \
    ARAGELI_ASSERT_0(p >= PREC_MIN && p <= PREC_MAX && "Precision is out of range")

#define CHECK_MODE(m)    \
    ARAGELI_ASSERT_0(m >= big_float_impl::exact_rounding && m <= big_float_impl::round_outward_zero && "Try to use incorrect rounding mode")

#if !defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE) ||    \
    defined(ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_BIG_FLOAT)

#include <cmath>
#include "big_float_impl.hpp"


namespace Arageli
{

template <typename T>
void big_float_impl::from_native_float(const T &f)
{
    typedef std::numeric_limits<T> Nl;

    //using namespace _Internal;

    ARAGELI_ASSERT_1(Nl::is_specialized);
    ARAGELI_ASSERT_1(!Nl::is_integer);
    ARAGELI_ASSERT_0(!(Nl::has_infinity && f == Nl::infinity()));
    ARAGELI_ASSERT_0(!(Nl::has_quiet_NaN && f == Nl::quiet_NaN()));
    ARAGELI_ASSERT_0(!(Nl::has_signaling_NaN && f == Nl::signaling_NaN()));
    ARAGELI_ASSERT_0(!(Nl::has_denorm && f == Nl::denorm_min()));

    if ( !f )
    {
        big_float_impl ();
        return;
    }

    int digits_need = (Nl::digits - 1) / _Internal::bits_per_digit;
    _Internal::digit *man;
    _Internal::digit *p = (_Internal::digit *)&f;

    man = new _Internal::digit [ digits_need + 1 ];
    ARAGELI_ASSERT_0(man && "the heap overflow");

    for ( int i = 0; i < digits_need; i++,p++ )
        man [ i ] = *p;

    int bits_remain = (Nl::digits - 1) % _Internal::bits_per_digit;

    man [ digits_need ] =
        (
            (*p << (_Internal::bits_per_digit - bits_remain)) >>
            (_Internal::bits_per_digit - bits_remain)
        ) | (1U << bits_remain);

    prec = Nl::digits;
    mode = big_float_impl::get_default_rounding_mode();

    #ifdef ARAGELI_GMP
    ARAGELI_ASSERT_ALWAYS_EX1(0, "GMP and big_float_impl aren't compatible");
    #else
    s.free_mem_and_alloc_number( f > 0 ? 1 : -1, man, digits_need + 1 );
    #endif

    int ex;
    frexp (f, &ex);
    e = ex - Nl::digits;
}

template <typename T>
void big_float_impl::from_native_int (const T& i)
{
    s = big_int(i);
    prec =
        (big_float_impl::get_default_precision() < std::numeric_limits<T>::digits) ?
        std::numeric_limits<T>::digits :
        big_float_impl::get_default_precision();

    CHECK_PREC(prec)
    mode = big_float_impl::get_default_rounding_mode();
}

template <typename T>
T big_float_impl::to_native_int () const
{
    return
        s.sign() > 0 ?
        ifloor(*this).to_native_int<int>() :
        iceil (*this).to_native_int<int>();
}

template <typename T>
T big_float_impl::to_native_float () const
{
    typedef std::numeric_limits<T> Nl;

    ARAGELI_ASSERT_1(Nl::is_specialized);
    ARAGELI_ASSERT_1(!Nl::is_integer);
    ARAGELI_ASSERT_0(Nl::has_infinity);

    big_int temp(s);
    int shift = temp.length() - Nl::digits;
    if (shift > 0 )
        temp >>= shift;

    /*!! */
    int tlen = temp.length();
    int digits_need = tlen / _Internal::bits_per_digit + 1;
    _Internal::digit *p = new _Internal::digit [ digits_need ];
    ARAGELI_ASSERT_0(p && "the heap overflow");

    #ifdef ARAGELI_GMP
    ARAGELI_ASSERT_ALWAYS_EX1(0, "GMP and big_float_impl aren't compatible");
    #else
    std::memmove( p, temp.number->data, digits_need * sizeof (_Internal::digit));
    #endif

    p [digits_need - 1] /*&= ~(_Internal::digit(1) << tlen %_Internal::bits_per_digit - 1)*/;
    T ret = s.sign() * ldexp (*((T*)p), e + (shift + Nl::max_exponent + Nl::digits - 3));
    delete [] p;
    return ret;
}

}//namespace Arageli

#else

#include <sstream>

#include "powerest.hpp"
#include "logarithm.hpp"
#include "big_float_impl.hpp"
#include <cmath>


namespace Arageli
{


// here: conversion from long double to size_t,
// possible loss of data
#pragma warning(disable : 4244)


//a-la bigarith
void big_float_impl_warning ( const char *s )
{
    std::cerr << "Big float warning: " << s << "\n";
}

void big_float_impl_error ( const char *s )
{
    std::cerr << "Big float error: " << s << "\n";
}

void big_float_impl_fatal_error(const char *s)
{
    std::cerr << "Big float fatal error: \n " << s << "\n";
    exit(1);
}

/*
* Constructors
*/
//constructor from big_int
big_float_impl::big_float_impl (const big_int& i)
{
    s = big_int(i);
    prec = i.length();
    CHECK_PREC(prec)
    mode = get_default_rounding_mode();
}

big_float_impl & big_float_impl::operator = (const big_int &b)
{
    s = b;
    normalize_1 ( prec, mode );
    return *this;
}

// Normalization for big_float_impl number
void big_float_impl::normalize_1 ( big_float_impl::prec_t prec, big_float_impl::rounding_mode_t mode )
{
    this->prec = prec;
    this->mode = mode;
    if ( !s.sign() ) //significant s = 0, then exponenta = 0
    {
        e = 0;
        return;
    }

    int s_len = s.length() - prec;

    if ( mode == exact_rounding && s.length() <= PREC_MAX ||  (!s_len)  )
        return;// do nothing - mantissa has apropriate length
    if ( s_len < 0  )
    {
        s = s << -s_len;
        e = e + s_len;
        return;
    }//write zero at the end and return
    ARAGELI_ASSERT_0
    (
        mode!= exact_rounding &&
        "Impossibly to perfom operarion with exact_rounding rounding mode ( length too long )"
    );

    /*
    if ( mode == exact_rounding )
    {
        mode = round_to_nearest; // s.lenght > PREC_MAX
        big_float_impl_error ( "can't perfom operarion with exact_rounding rounding mode" );
    }
    */

    int is_digit = 0; // needs then rounds to +-infinity
    is_digit = ( (s >> s_len - 1) << s_len - 1 != s );
    s = s >> s_len - 1;
    e = e + s_len - 1;

    switch (mode)
    {
        case round_to_nearest:
                        s = ( s + s.sign() ) >> 1;
                        e = e + 1;
                        if ( s.length() > prec )
                        {
                            s = s >> 1;
                            e = e + 1;
                        }
                        break;
        case round_toward_zero:
                        s = s >> 1;
                        e = e + 1;
                        break;
        case round_toward_infinity:
                        if ( s.sign() < 0 )
                            s = s >> 1;
                        else
                        {
                            if ( is_digit )
                                s = ( s >> 1) + 1;
                            else
                                s = s + 1 >> 1;
                        }
                        e = e + 1;

                        if ( s.length() > prec )
                        {
                            s = s >> 1;
                            e = e + 1;
                        }
                        break;
        case round_toward_neg_infinity:
                        if ( s.sign() > 0  )
                            s = s >> 1;
                        else
                        {
                            if ( is_digit )
                                s = ( s >> 1) - 1;
                            else
                                s = s - 1 >> 1;
                        }
                        e = e + 1;

                        if ( s.length() > prec )
                        {
                            s = s >> 1;
                            e = e + 1;
                        }
                        break;
        case round_outward_zero:
                        //s=s/(max_digit+1);
                        if ( is_digit )
                            s = ( s >> 1 ) + s.sign();
                        else
                        {
                            if ( s [ 0 ] )
                                s = ( s >> 1 ) + s.sign();
                            else
                                s = s >> 1;
                        }
                        e = e + 1;

                        if ( s.length() > prec )
                        {
                            s = s >> 1;
                            e = e + 1;
                        }
                        break;
        default:
                        break;

    }//switch(mode)
    return;
}


/*
*  compare operations
*/

int cmp (const big_float_impl & a, const big_float_impl & b)
{
    big_float_impl temp (sub(a, b, std::max(a.get_precision(), b.get_precision()), big_float_impl::round_to_nearest));
    return temp.sign();
}

//adding with precision prec, rounding mode is mode
//(it is all the same the b and c are normlized or not)
//TODO remove WARNINGS related with exact_rounding mode
big_float_impl add ( const big_float_impl & b, const big_float_impl & c, big_float_impl::prec_t prec, big_float_impl::rounding_mode_t mode )
{
    CHECK_PREC(prec)
    CHECK_MODE(mode)
    big_int b_e(b.e), c_e(c.e), b_s(b.s), c_s(c.s);

    big_float_impl temp;

    // some checks
    if ( !c_s.length() ) // c == 0, so b + c == b
    {
        temp = b;
        temp.normalize_1 ( prec, mode );
        return temp;
    }
    if ( !b_s.length() ) // b == 0, so b + c == c
    {
        temp = c;
        temp.normalize_1 ( prec, mode );
        return temp;
    }

    // followed if ( ... )... needed for adding no normalized numbers, make equal length
    long delta = b_s.length() - c_s.length();

    if ( delta > 0 )
    {
        c_s = c_s << delta;
        c_e = c_e - delta;
    }
    else
    {
        b_s = b_s << -delta;
        b_e = b_e + delta;
    }

    if ( b_e < c_e )
    {
        swap ( b_e, c_e );
        swap ( b_s, c_s );
    }

    big_int ed = b_e - c_e;//or three times to compute or once to save

    if ( ed > prec  && mode != big_float_impl::exact_rounding)
    {
        temp.s = b_s;//
        temp.e = b_e;
        if ( mode != big_float_impl::round_to_nearest )
        {
            temp.s = ( temp.s << 1 ) + c_s.sign();
            temp.e = temp.e - 1;
        }
    }
    else
    {
        if ( ed.length() > _Internal::bits_per_digit /* PREC_MAX - b_s.length()*/ ) //ed > PREC_MAX and mode == exact_rounding
            big_float_impl_fatal_error( "can't perfom adding with exact_rounding rounding mode" );

        b_s = b_s << ed/*.to_digit ()*/;

        temp.s = c_s + b_s;
        temp.e = c_e;
    }
    temp.normalize_1 ( prec, mode );
    return temp;
}

//subtraction
big_float_impl sub(const big_float_impl & b, const big_float_impl & c, big_float_impl::prec_t prec, big_float_impl::rounding_mode_t mode)
{
    big_float_impl temp(c);
    temp.s = -temp.s;
    temp = add (b, temp, prec, mode);
    return temp;
}

//myltiplying with precision prec, rounding mode is mode
//(it is all the same are the b and c normlized or not)
big_float_impl mul(const big_float_impl & b, const big_float_impl & c, big_float_impl::prec_t prec, big_float_impl::rounding_mode_t mode)
{
    CHECK_PREC(prec)
    CHECK_MODE(mode)
    big_float_impl temp;
    //TODO The way to simplify computation is check prec and then truncate b or/and c
    temp.s = b.s * c.s;
    temp.e = b.e + c.e;
    temp.normalize_1 ( prec, mode );
    return temp;
}

//division (temporary version)
big_float_impl div (const big_float_impl & b, const big_float_impl & c, big_float_impl::prec_t prec, big_float_impl::rounding_mode_t mode)
{
    if ( c.s==0 )
        return big_float_impl();//
    big_float_impl x;
    big_int temp = pow2<big_int> ( prec << 1 );

    x.s = temp / c.s;
    x.e = -c.e - prec - prec;

    x.normalize_1 (prec,mode );
    return mul( x, b, prec, mode );
}
//not used for the time present
//big_float_impl divnu (const big_float_impl & b, const big_float_impl & c, big_float_impl::prec_t prec, big_float_impl::rounding_mode_t mode)
//{
//    short k = 0;
//    big_float_impl x;
//    big_float_impl temp = 1;
//    temp.out ( std::cout, 'd');
//    x.s = c.s;
//    while ( x.s.length() % _Internal::bits_per_digit > 0 )
//    {
//        x.s = x.s << 1;
//        k++;
//    }
//    x.e = -( int )x.s.length();// / bits_per_digit;
//    for ( int i=0; i < prec + 1; i++ )
//    {
//        temp = mul (2,temp,prec,mode) - mul ( mul ( temp,temp,prec,mode ),x,prec,mode );
//    }
//    for (int i = 0; i < k; i++)
//    {
//        temp = temp * big_float_impl(2.0);
//    }
//    temp.e = temp.e + x.e - c.e;
//    temp = mul ( b, temp, prec, mode );
//    return temp;
//}

//experimental method
big_float_impl div_i ( const big_int & c )
{
    std::size_t l  = c.length();

    std::size_t step = 1;
    big_float_impl res;

    res.s = 1;

    for ( step = 1 ; step < 32; step <<= 1 )
    {
        res.s = res.s * ( pow2<big_int> ( l * step + 1  ) - c * res.s );
    }
    //res.e =- l * ( step );
    std::cout << 10000 * res.s / pow2<big_int> ( 500)<< std::endl;
    //std::cout << res.s << std::endl << res.e << std::endl;
    return res;
}

std::pair<std::basic_string<char>, big_float_impl::exp_t> big_float_impl::to_string (std::size_t n, int /*base*/) const
{
    big_int dm, de;
    do_dec_convert(dm, de, n - 1, s,e);
    std::ostringstream str;
    str << dm;
    return std::make_pair(str.str(), de + str.str().length() - (str.str()[0] == '-'));
}

/*
* miscellaneous functions
*/
//returns a floating-point value representing the smallest integer that is greater than or equal to a
big_float_impl ceil  (const big_float_impl & a)
{
    big_float_impl res ( a );
    big_int e = a.e;
    big_int temp;

    if ( e >=  0 )
        return res;
    if ( e <= -big_int(a.s.length()) )
    {
        if ( a.s > 0 )
            return 1.0;
        else
            return 0.0;
    }
    temp = ( res.s >> -e) << -e;
    if ( res.s != temp && res.s > 0 )
        res.s = temp + pow2<big_int>( -e/*.to_digit ()*/ );
    else
        res.s = temp;
    return res;
}

// returns a floating-point value representing the largest integer that is less than or equal to a
big_float_impl floor (const big_float_impl & a)
{
    big_float_impl res ( a );
    big_int e = a.e;
    big_int temp;

    if ( e >=  0 )
        return res;
    if ( e <= -big_int(a.s.length()) )
    {
        if ( a.s > 0 )
            return 1.0;
        else
            return 0.0;
    }
    e = -e;
    temp = ( res.s >> e) << e;
    if ( res.s != temp && res.s < 0 )
        res.s = temp - pow2<big_int>( e );
    else
        res.s = temp;
    return res;

}

// Returns the fractal part of the number
big_float_impl frac  (const big_float_impl & a)
{
    if ( a.e >= 0 )
        return 0.0;

    big_float_impl res ( a );
    if ( res.e <= -big_int(res.s.length()) )
        return res;

    res.s = res.s - ((res.s >> -res.e) << -res.e);
    return res;
}

// returns an big-integer value representing the smallest integer that is greater than or equal to a
big_int iceil (const big_float_impl & a)
{
    std::size_t l = a.s.length();
    if ( !l )
        return l;

    big_int temp ( a.s );
    if ( a.e > 0 )
    {
        //if ( a.e + a.s.length() > SHRT_MAX/*??*/  )
        //    big_float_impl_fatal_error ( "Can't convert from big_float_impl to big_int ( exponent too large )..." );
        return a.s <<a.e;
    }
    if ( a.e <= -big_int(a.s.length()) )
        return ( a.s > 0 ) ? 1 : 0;
    temp = temp >> -a.e;
    if ( temp << -a.e != a.s && temp > 0)
        return temp + 1;
    return temp;
}

// returns a big-integer  value representing the largest integer that is less than or equal to a
big_int ifloor(const big_float_impl & a)
{
    std::size_t l = a.s.length();
    if ( !l )
        return l;

    big_int temp ( a.s );
    if ( a.e > 0 )
    {
        //if ( a.e + a.s.length() > SHRT_MAX/*??*/  )
        //    big_float_impl_fatal_error ( "Can't convert from big_float_impl to big_int ( exponent too large )..." );
        return a.s <<a.e;
    }
    if ( a.e <= -big_int(a.s.length()) )
        return ( a.s > 0 ) ? 0 : -1;
    temp = temp >> -a.e;
    if ( temp << -a.e != a.s && temp > 0)
        return temp;
    return temp - 1;
}

big_float_impl big_float_impl::sqrt (big_float_impl::prec_t p, big_float_impl::rounding_mode_t m) const
{
    if (s == 0)
        return big_float_impl (p, m);
    big_float_impl res (*this);
    int l;

    if (res.e [0])
    {
        res.s = res.s << 1;
        res.e = res.e - 1;
    }

    l = (res.s.length() - 2 * prec - 1);
    l = (l / 2) * 2;

    if (l > 0)
    {
        res.s = res.s >> l;
        res.e = res.e + l;
    }
    else
    {
        res.s = res.s << (-l) + 2;
        res.e = res.e + l - 2;
    }
    res.s = Arageli::sqrt (res.s);
    res.e = res.e >> 1;
    res.normalize_1 (prec,mode);//error may be if mode == round_outward_zero

    return res;
}

//Newton fsqrt
big_float_impl nfsqrt ( const big_float_impl & bf, big_float_impl::prec_t prec, big_float_impl::rounding_mode_t mode )
{
    big_float_impl res;
    //std::cout << bf << std::endl;
    res.e = ((( bf.get_exp() ) + bf.get_significant().length() - 1) >> 1) + 1;
    res.s = 1;
    //res.normalize_1();
//    std::cout << "x0 "<<res << std::endl;
    //res.out ( std::cout, 'd' ) ;
    //std::cout << std::endl;
    std::size_t n = (long) (log ( (long double)prec + 1 ) / log ( 2.0l ) + 1.0);

    for ( std::size_t counter = 1; counter <= 2 * n; counter ++ )
    {
//        std::cout << div( bf, res, counter * 2, big_float_impl::rounding_mode_t::exact_rounding ) << std::endl;
        res = add ( res,  div( bf, res, counter * 2, big_float_impl::exact_rounding ), counter * 2, big_float_impl::round_to_nearest );
//        std::cout << res << std::endl;
        res.e = res.e - 1;
//       std::cout << counter << " iter \t";
//       std::cout << res << std::endl;
        //res.out( std::cout, 'd' );
        //std::cout << std::endl;
    }
    res.normalize_1 (prec, mode);//error may be if mode == round_outward_zero

    return res;
}
//Random numbers
big_float_impl frandom  ( long prec )
{
    CHECK_PREC(prec)//
    return big_float_impl (big_int::random_with_length_or_less(prec), big_int(-prec));
}
//Random numbers
big_float_impl frandom  ( )
{
    return big_float_impl ( big_int::random_with_length_or_less(big_float_impl::get_default_precision() ), big_int (-big_float_impl::get_default_precision()) );
}

big_float_impl frandom1 ( long bits, const big_int & exp)
{
    CHECK_PREC(bits)//
    return big_float_impl (  big_int::random_with_length_or_less( bits ) + pow2<big_int>( bits - 1 ) , exp );
}

}//namespace Arageli
#endif

