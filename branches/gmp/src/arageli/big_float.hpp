/*****************************************************************************
*
* big_float.hpp
*
* This file is a part of the Arageli library.
*
* Copyright (C) 2008 Alexander Pshenichnikov
* Copyright (C) 2008 Nikolay Santalov
* University of Nizhni Novgorod, Russia
*
* The Arageli Library is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License version 2
* as published by the Free Software Foundation.
*
* The Arageli Library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*
* We are also open for dual licensing for the whole library or
* for its particular part. If you are interested to get the library
* in this way, i.e. not under the GNU General Public License,
* please contact Arageli Support Service support.arageli@gmail.com.
*
******************************************************************************/

/**
*     \file big_float.hpp
*     \brief Big float number class implementation.
*
*     This module implements a class big_float_t for representing
*     Big Float Numbers, i.e. unlimited precision real numbers with
*     floating point.
*
**/

//TODO really need to include all these headers here? Think about moving corresponding functionality to separate files
#include <cmath> //for isnan function.

#include <iostream> //for ostream and istream
#include <iomanip> //for setw, setfill, etc.
#include <sstream> //for stringstream
#include <utility> //for pair
#include <string> //for basic_string

#include "config.hpp"

namespace Arageli
{

template <typename BigFloatImpl>
class big_float_t
{
    typedef BigFloatImpl bf_t;

public:
    typedef bf_t rep_t;
    typedef typename rep_t::prec_t prec_t;
    typedef typename rep_t::exp_t exp_t;

    typedef typename rep_t::rounding_mode_t rounding_mode_t;

    //constructors
    big_float_t ()
    {}

    big_float_t (prec_t prec, rounding_mode_t m) :
    rep (prec, m)
    {}

    template <typename T>
    big_float_t (T x) :
    rep (x)
    {}

    big_float_t (const big_float_t<rep_t> &b) :
    rep(b.get_rep())
    {}

    ~big_float_t ()
    {}

    //operator =
    big_float_t & operator = (const big_float_t<rep_t> &other)
    {
        if (this != &other)
        {
            rep = other.get_rep();
        }
        return *this;
    }

    template<typename T>
    big_float_t & operator = (T x)
    {
        rep = x;
        return *this;
    }

    //conversion
    template <typename T>
    operator T () const
    {
        return T(get_rep());
    }

    //precision and round mode
    static prec_t get_default_precision()
    {
        return rep_t::get_default_precision();
    }

    static rounding_mode_t get_default_rounding_mode()
    {
        return rep_t::get_default_rounding_mode();
    }

    static void set_default_precision()
    {
        return rep_t::set_default_precision();
    }

    static void set_default_rounding_mode()
    {
        return rep_t::set_default_rounding_mode();
    }

    void set_precision(prec_t prec)
    {
        return get_rep().set_precision(prec);
    }

    prec_t get_precision() const
    {
        return get_rep().get_precision();
    }

    rounding_mode_t get_rounding_mode () const
    {
        return get_rep().get_rounding_mode();
    }

    void set_rounding_mode (rounding_mode_t m)
    {
        return get_rep().set_rounding_mode(m);
    }

    //representation
    const rep_t & get_rep () const
    {
        return rep;
    }

    rep_t & get_rep ()
    {
        return rep;
    }

    exp_t get_exp() const
    {
        return get_rep().get_exp();
    }

    void setsign (int sign)
    {
        get_rep().setsign(sign);
    }

private:
    rep_t rep;
};

template <typename T>
inline bool isnan(const T)
{
    return false;
}

template <>
inline bool isnan(const float f)
{
    #ifdef _MSC_VER
    return static_cast<bool>(::_isnan(f));
    #else
        return static_cast<bool> (std::isnan(f));
    #endif
}

template <>
inline bool isnan(const double f)
{
    #ifdef _MSC_VER
        return static_cast<bool>(::_isnan(f));
    #else
        return static_cast<bool> (std::isnan(f));
    #endif
}

template <>
inline bool isnan(const long double f)
{
    #ifdef _MSC_VER
        return static_cast<bool>(::_isnan(f));
    #else
        return static_cast<bool> (std::isnan(f));
    #endif
}

//comparison
template <typename BF, typename T>
inline int cmp (const big_float_t<BF> &b, const T c)
{
    return cmp(b.get_rep(), c);
}

template <typename BF>
inline int cmp (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return cmp(b.get_rep(), c.get_rep());
}

template <typename BF, typename T>
inline int cmp (const T &c, const big_float_t<BF> &b)
{
    return -cmp(b.get_rep(), c);
}

template <typename BF, typename T>
inline bool operator == (const big_float_t<BF> &b, const T x)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) == 0;
}

template <typename BF, typename T>
inline bool operator > (const big_float_t<BF> &b, const T x)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) > 0;
}

template <typename BF, typename T>
inline bool operator < (const big_float_t<BF> &b, const T x)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) < 0;
}

template <typename BF, typename T>
inline bool operator <= (const big_float_t<BF> &b, const T x)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) <= 0;
}

template <typename BF, typename T>
inline bool operator >= (const big_float_t<BF> &b, const T x)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) >= 0;
}

template <typename BF, typename T>
inline bool operator != (const big_float_t<BF> &b, const T x)
{
    return  (isnan(b) || isnan(x)) ? true : cmp(b,x) != 0;
}

template <typename BF, typename T>
inline bool operator == (const T x, const big_float_t<BF> &b)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) == 0;
}

template <typename BF, typename T>
inline bool operator > (const T x, const big_float_t<BF> &b)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) < 0;
}

template <typename BF, typename T>
inline bool operator < (const T x, const big_float_t<BF> &b)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) > 0;
}

template <typename BF, typename T>
inline bool operator <= (const T x, const big_float_t<BF> &b)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) >= 0;
}

template <typename BF, typename T>
inline bool operator >= (const T x, const big_float_t<BF> &b)
{
    return  (isnan(b) || isnan(x)) ? false : cmp(b,x) <= 0;
}

template <typename BF, typename T>
inline bool operator != (const T x, const big_float_t<BF> &b)
{
    return  (isnan(b) || isnan(x)) ? true : cmp(b,x) != 0;
}

template <typename BF>
inline bool operator == (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return  (isnan(b) || isnan(c)) ? false : cmp(b,c) == 0;
}

template <typename BF>
inline bool operator > (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return  (isnan(b) || isnan(c)) ? false : cmp(b,c) > 0;
}

template <typename BF>
inline bool operator < (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return  (isnan(b) || isnan(c)) ? false : cmp(b,c) < 0;
}

template <typename BF>
inline bool operator >= (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return  (isnan(b) || isnan(c)) ? false : cmp(b,c) >= 0;
}

template <typename BF>
inline bool operator <= (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return  (isnan(b) || isnan(c)) ? false : cmp(b,c) <= 0;
}

template <typename BF>
inline bool operator != (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return  (isnan(b) || isnan(c)) ? true : cmp(b,c) != 0;
}

template <typename BF>
inline big_float_t<BF> operator + (const big_float_t<BF> &b)
{
    return b;
}

template <typename BF>
inline big_float_t<BF> operator - (const big_float_t<BF> &b)
{
    big_float_t<BF> ret(b);
    ret.setsign(-1);
    return ret;
}

//arithmetic operations (+, *, -, /)
template <typename BF, typename T>
inline big_float_t<BF> add
(
    const big_float_t<BF> &b,
    T x,
    typename big_float_t<BF>::prec_t prec = big_float_t<BF>::get_default_precision,
    typename big_float_t<BF>::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(add(b.get_rep(), x, prec, mode));
}

template <typename BF, typename T>
inline big_float_t<BF> add
(
    T x,
    const big_float_t<BF> &b,
    typename big_float_t<BF>::prec_t prec = big_float_t<BF>::get_default_precision,
    typename big_float_t<BF>::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(add(b.get_rep(), x, prec, mode));
}

template <typename BF>
inline big_float_t<BF> add
(
    const big_float_t<BF> &b,
    const big_float_t<BF> &c,
    typename big_float_t<BF>::prec_t prec = big_float_t<BF>::get_default_precision,
    typename big_float_t<BF>::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(add(b.get_rep(), c.get_rep(), prec, mode));
}

template <typename BF, typename T>
inline big_float_t<BF> mul
(
    const big_float_t<BF> &b,
    T x,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(mul(b.get_rep(), x, prec, mode));
}

template <typename BF, typename T>
inline big_float_t<BF> mul
(
    T x,
    const big_float_t<BF> &b,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(mul(b.get_rep(), x, prec, mode));
}

template <typename BF>
inline big_float_t<BF> mul
(
    const big_float_t<BF> &b,
    const big_float_t<BF> &c,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(mul(b.get_rep(), c.get_rep(), prec, mode));
}

template <typename BF, typename T>
inline big_float_t<BF> sub
(
    const big_float_t<BF> &b,
    T x,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(sub(b.get_rep(), x, prec, mode));
}

template <typename BF, typename T>
inline big_float_t<BF> sub
(
    T x,
    const big_float_t<BF> &b,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(sub(x, b.get_rep(), prec, mode));
}

template <typename BF>
inline big_float_t<BF> sub
(
    const big_float_t<BF> &b,
    const big_float_t<BF> &c,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(sub(b.get_rep(), c.get_rep(), prec, mode));
}

template <typename BF, typename T>
inline big_float_t<BF> div
(
    const big_float_t<BF> &b,
    T x,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(div(b.get_rep(), x, prec, mode));
}

template <typename BF, typename T>
inline big_float_t<BF> div
(
    T x,
    const big_float_t<BF> &b,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(div(x, b.get_rep(), prec, mode));
}

template <typename BF>
inline big_float_t<BF> div
(
    const big_float_t<BF> &b,
    const big_float_t<BF> &c,
    typename BF::prec_t prec = big_float_t<BF>::get_default_precision(),
    typename BF::rounding_mode_t mode = big_float_t<BF>::get_default_rounding_mode()
)
{
    return big_float_t<BF>(div(b.get_rep(), c.get_rep(), prec, mode));
}

//operator +
template<typename BF, typename T>
inline big_float_t<BF> operator + (const big_float_t<BF> &b, T x)
{
    return add(b, x, b.get_precision());
}

template<typename BF, typename T>
inline big_float_t<BF> operator + (T x, const big_float_t<BF> &b)
{
    return add(b, x, b.get_precision());
}

template<typename BF>
inline big_float_t<BF> operator + (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return add(b, c, std::max(b.get_precision(), c.get_precision()));
}

//operator *
template<typename BF, typename T>
inline big_float_t<BF> operator * (const big_float_t<BF> &b, T x)
{
    return mul(b, x, b.get_precision());
}

template<typename BF, typename T>
inline big_float_t<BF> operator * (T x, const big_float_t<BF> &b)
{
    return mul(b, x, b.get_precision());
}

template<typename BF>
inline big_float_t<BF> operator * (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return mul(b, c, std::max(b.get_precision(), c.get_precision()));
}

//operator -
template<typename BF, typename T>
inline big_float_t<BF> operator - (const big_float_t<BF> &b, T x)
{
    return sub(b, x, b.get_precision());
}

template<typename BF, typename T>
inline big_float_t<BF> operator - (T x, const big_float_t<BF> &b)
{
    return sub(x, b, b.get_precision());
}

template<typename BF>
inline big_float_t<BF> operator - (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return sub(b, c, std::max(b.get_precision(), c.get_precision()));
}

//operator /
template<typename BF, typename T>
inline big_float_t<BF> operator / (const big_float_t<BF> &b, T x)
{
    return div(b, x, b.get_precision());
}

template<typename BF, typename T>
inline big_float_t<BF> operator / (T x, const big_float_t<BF> &b)
{
    return div(x, b, b.get_precision());
}

template<typename BF>
inline big_float_t<BF> operator / (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return div(b, c, std::max(b.get_precision(), c.get_precision()));
}

//special numbers
template <typename BF>
inline bool issingular (const big_float_t<BF> &b)
{
    return b.get_rep().issingular();
}

template <typename BF>
inline bool isnan (const big_float_t<BF> &b)
{
    return b.get_rep().isnan();
}

template <typename BF>
inline bool isinf (const big_float_t<BF> &b)
{
    return b.get_rep().isinf();
}

template <typename BF>
inline bool isminf (const big_float_t<BF> &b)
{
    return b.get_rep().isminf();
}

template <typename BF>
inline bool ispinf (const big_float_t<BF> &b)
{
    return b.get_rep().ispinf();
}

template <typename BF>
inline bool isfinite (const big_float_t<BF> &b)
{
    return b.get_rep().isfinite();
}

template <typename BF>
inline bool iszero (const big_float_t<BF> &b)
{
    return b.get_rep().iszero();
}

template <typename BF>
inline bool ismzero (const big_float_t<BF> &b)
{
    return b.get_rep().ismzero();
}

template <typename BF>
inline bool ispzero (const big_float_t<BF> &b)
{
    return b.get_rep().ispzero();
}

//TEMPORARY VERSION!
template <typename BF, typename T>
inline big_float_t<BF> operator += (big_float_t<BF> &b, T x)
{
    return (b = b + x);
}

template <typename BF, typename T>
inline big_float_t<BF> operator -= (big_float_t<BF> &b, T x)
{
    return (b = b - x);
}

template <typename BF, typename T>
inline big_float_t<BF> operator /= (big_float_t<BF> &b, T x)
{
    return (b = b / x);
}

template <typename BF, typename T>
inline big_float_t<BF> operator *= (big_float_t<BF> &b, T x)
{
    return (b = b * x);
}

template <typename BF>
inline big_float_t<BF> operator -= (big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return (b = b + c);
}

template <typename BF>
inline big_float_t<BF> operator /= (big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return (b = b / c);
}

template <typename BF>
inline big_float_t<BF> operator *= (big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return (b = b * c);
}

//special functions

//miscelanious
template <typename BF>
inline int sign(const big_float_t<BF> &b)
{
    return b.get_rep().sign();
}

//output and input
namespace _Internal
{

template <typename BF, typename Ostream>
Ostream & big_float_t_out_scientific (Ostream &os, std::ios_base::fmtflags flags, const big_float_t<BF> &b)
{
    std::pair <std::basic_string<char>, typename BF::exp_t> bs = b.get_rep().to_string(os.precision(), 10);
    std::basic_string<char> &man = bs.first;
    typename BF::exp_t &exp = bs.second;
    os << man.insert(1 + (sign(b) < 0), ".") << (flags & std::ios_base::uppercase ? 'E' : 'e') << exp - 1;
    return os;
}

template<typename BF, typename Ostream>
Ostream & big_float_t_out_auto (Ostream &os, std::ios_base::fmtflags flags, const big_float_t<BF> &b)
{
    using namespace std;
    struct zero_eraser
    {
        std::basic_string<char> & operator () (std::basic_string<char> &str)
        {
            if (str.find('.') != std::basic_string<char>::npos)
            {
                std::size_t zpos = str.find_last_not_of('0');
                str.erase (str[zpos] == '.' ? zpos: zpos + 1);
            }
            return str;
        }
    };

    long precision = os.precision();
    pair <basic_string<char>, typename BF::exp_t> bs = b.get_rep().to_string(precision, 10);
    basic_string<char> &man = bs.first;

    if (man[0] == '0')
    {
        os << "0.0";
        return os;
    }
    else
    {
        typename BF::exp_t &exp = bs.second;
        int sigm = sign(b) < 0;
        zero_eraser ze;

        if (exp > precision || exp < -3)
        {
            os << ze(man.insert(1 + sigm, ".")) << (flags & ios_base::uppercase ? 'E' : 'e') << exp - 1;
            return os;
        }
        else if ( exp <= 0 )
        {
            stringstream ss;
            os << ze(man.insert(sigm, static_cast<stringstream &>(ss << setfill('0') << setw(-exp+2) << left << "0.").str()));
            return os;
        }
        else if ( exp < precision )
        {
            os << ze(man.insert(exp + sigm , "."));
            return os;
        }
        os << man;
        return os;
    }
}

#if !defined (M_LOG10E)
    #define M_LOG10E 0.434294481903251827651
#endif

#if !defined (M_LOG2E)
    #define M_LOG2E 0.434294481903251827651
#endif

#define M_LOG10_2 M_LOG10E / M_LOG2E

template <typename BF, typename Ostream>
Ostream & big_float_t_out_fixed ( Ostream &os, std::ios_base::fmtflags flags, const big_float_t<BF> &b)
{
    using namespace std;
    typename BF::exp_t exp = b.get_exp();
    //approximate how many decimal digits we need to retrieve
    size_t precision = os.precision();
#if(__STDC_VERSION__ >= 199901L)
    long position = 1 + precision + lround (static_cast<double> (exp) * M_LOG10_2);
#else
    long position = 1 + precision + ceil (static_cast<double> (exp) * M_LOG10_2);
#endif
    if ( position <= 0)
    {
        os << "0.0";
        return os;
    }

    if (position == 1) ++position;

    pair <basic_string<char>, typename BF::exp_t> bs = b.get_rep().to_string(position, 10);
    basic_string<char> &man = bs.first;
    exp = bs.second;

    int sigm = sign(b) < 0;

    stringstream ss;
    if (exp <= 0)
        man.insert(sigm, static_cast<stringstream &> (ss << "0." << setw(-exp) << setfill('0') << "").str());
    else
        man.insert(exp + sigm, ".");
    man.erase (man.end()-(position - exp - precision), man.end());
    os << man;
    return os;
}

}//namespace _Internal

template <typename BF, typename Ostream>
Ostream & operator << (Ostream &os, const Arageli::big_float_t<BF> &b)
{
    if ((os.flags() & std::ios_base::showpos) && sign(b) > 0) os << '+';

    std::ios_base::fmtflags fl = os.flags();

    if (fl & std::ios_base::scientific)
    {
        return _Internal::big_float_t_out_scientific (os, fl, b);
    }
    else if ( fl & std::ios_base::fixed)
    {
        return _Internal::big_float_t_out_fixed (os, fl, b);
    }

    else
    {
        return _Internal::big_float_t_out_auto (os, fl, b);
    }

    return os;
}

template <typename BF, typename Stream>
Stream & operator >> (Stream &is, big_float_t<BF> &b)
{
    ARAGELI_ASSERT_ALWAYS("std::istream & operator >> (std::istream &is, big_float_t<BF> &b) not implemenetd yet!")
    //WARNING! temporary solution! TODO add input error handling
    //std::basic_string<char> str;
    //is >> str;
    //mpfr_strtofr(b.get_rep(), str.c_str(), 0, 10, b.get_mp_rounding_mode());
    return is;
}

//special functions
#define _ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(FUNC)    \
template <typename BF>    \
inline big_float_t<BF> FUNC (const big_float_t<BF> &b, typename big_float_t<BF>::prec_t prec, typename big_float_t<BF>::rounding_mode_t mode)    \
{    \
    return big_float_t<BF>(/*FUNC(b.get_rep(), prec, mode)*/);    \
}    \

#define _ARAGELI_BIG_FLOAT_SPECIAL_BFUNC(FUNC)    \
template <typename BF>    \
inline big_float_t<BF> FUNC    \
(    \
    const big_float_t<BF> &b,    \
    const big_float_t<BF> &c,    \
    typename big_float_t<BF>::prec_t prec,    \
    typename big_float_t<BF>::rounding_mode_t mode    \
)    \
{    \
    return big_float_t<BF>(FUNC(b.get_rep(), c.get_rep(), prec, mode));    \
}    \

#define _ARAGELI_BIG_FLOAT_SPECIAL_MIXED_BFUNC(FUNC)    \
template <typename BF>    \
inline big_float_t<BF> FUNC    \
(    \
    long x,    \
    const big_float_t<BF> &b,    \
    typename big_float_t<BF>::prec_t prec,    \
    typename big_float_t<BF>::rounding_mode_t mode    \
)    \
{    \
    return big_float_t<BF>(FUNC(x, b.get_rep(), prec, mode));    \
}    \

#define _ARAGELI_BIG_FLOAT_SPECIAL_VFUNC(FUNC)    \
template <typename BF>    \
inline big_float_t<BF> FUNC (typename big_float_t<BF>::prec_t prec, typename big_float_t<BF>::rounding_mode_t mode)    \
{    \
    return big_float_t<BF>(FUNC(prec, mode));    \
}    \

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(log)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(log2)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(log10)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(exp)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(exp2)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(exp10)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(sin)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(cos)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(tan)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(sec)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(csc)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(cot)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(asin)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(acos)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(atan)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(sinh)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(cosh)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(tanh)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(sech)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(csch)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(coth)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(asinh)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(acosh)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(atanh)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(log1p)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(expm1)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(eint)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(erf)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(erfc)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(gamma)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(lngamma)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(zeta)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(j0)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(j1)
_ARAGELI_BIG_FLOAT_SPECIAL_MIXED_BFUNC(jn)

_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(y0)
_ARAGELI_BIG_FLOAT_SPECIAL_UFUNC(y1)
_ARAGELI_BIG_FLOAT_SPECIAL_MIXED_BFUNC(yn)

_ARAGELI_BIG_FLOAT_SPECIAL_VFUNC(pi)
_ARAGELI_BIG_FLOAT_SPECIAL_VFUNC(log2)
_ARAGELI_BIG_FLOAT_SPECIAL_VFUNC(euler)
_ARAGELI_BIG_FLOAT_SPECIAL_VFUNC(catalan)

_ARAGELI_BIG_FLOAT_SPECIAL_BFUNC(atan2)
_ARAGELI_BIG_FLOAT_SPECIAL_BFUNC(hypot)

#undef _ARAGELI_BIG_FLOAT_SPECIAL_UFUNC
#undef _ARAGELI_BIG_FLOAT_SPECIAL_BFUNC
#undef _ARAGELI_BIG_FLOAT_SPECIAL_VFUNC
#undef _ARAGELI_BIG_FLOAT_SPECIAL_MIXED_BFUNC

template<typename BF>
inline big_float_t<BF> fac (unsigned long x, typename big_float_t<BF>::prec_t prec, typename big_float_t<BF>::rounding_mode_t mode)
{
    return big_float_t<BF> (fac(x, prec, mode));
}

template<typename BF>
inline big_float_t<BF> zeta (unsigned long x, typename big_float_t<BF>::prec_t prec, typename big_float_t<BF>::rounding_mode_t mode)
{
    return big_float_t<BF> (zeta, prec, mode);
}

//miscellaneous
template <typename BF>
void swap (big_float_t<BF> &b, big_float_t<BF> &c)
{
    return b.get_rep().swap(c.get_rep());
}

template<typename BF>
int signbit(const big_float_t<BF> &b)
{
    return b.signbit();
}

template <typename BF>
big_float_t<BF> abs(const big_float_t<BF> &b)
{
    return b.get_rep().abs();
}

//integer related functions
#define _ARAGELI_BIG_FLOAT_MISC_FUNC(FUNC)    \
template <typename BF>    \
inline big_float_t<BF> FUNC (const big_float_t<BF> &b)    \
{    \
    return big_float_t<BF>(FUNC(b));    \
}    \

_ARAGELI_BIG_FLOAT_MISC_FUNC(rint)
_ARAGELI_BIG_FLOAT_MISC_FUNC(ceil)
_ARAGELI_BIG_FLOAT_MISC_FUNC(floor)
_ARAGELI_BIG_FLOAT_MISC_FUNC(round)
_ARAGELI_BIG_FLOAT_MISC_FUNC(trunc)

_ARAGELI_BIG_FLOAT_MISC_FUNC(frac)
_ARAGELI_BIG_FLOAT_MISC_FUNC(nextbelow)
_ARAGELI_BIG_FLOAT_MISC_FUNC(nextabove)

#undef _ARAGELI_BIG_FLOAT_MISC_FUNC

template <typename BF>
inline big_float_t<BF> nexttoward (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    return big_float_t<BF>(nexttoward(b, c));
}


//min and max
template <typename BF>
big_float_t<BF> max (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    if (issingular(b) || issingular(c))
    {
        if (isnan(b) && isnan(c))
            return b;
        else if (isnan(b))
            return c;
        else if (isnan(c))
            return b;
        else if (iszero(b) && iszero(c))
        {
            if (ismzero(b))
                return c;
            else
                return b;
        }
    }
    if (cmp(b,c) <= 0)
        return c;
    else
        return b;
}

template <typename BF>
big_float_t<BF> min (const big_float_t<BF> &b, const big_float_t<BF> &c)
{
    if (issingular(b) || issingular(c))
    {
        if (isnan(b) && isnan(c))
            return b;
        else if (isnan(b))
            return c;
        else if (isnan(c))
            return b;
        else if (iszero(b) && iszero(c))
        {
            if (ismzero(b))
                return c;
            else
                return b;
        }
    }
    if (cmp(b,c) <= 0)
        return b;
    else
        return c;
}

template <typename Outiter, typename BF>
inline void generate_range_helper (big_float_t<BF>& t1, const big_float_t<BF>& t2, Outiter outiter)
{
    generate_range_helper_wo_inc(t1, t2, outiter);
}

}//namespace Arageli

#if defined(ARAGELI_MPFR)
    #include "mpfr_wrapper.hpp"
    namespace Arageli
    {
        typedef Arageli::big_float_t<Arageli::mpfr_wrapper> big_float;
    }
#else
    #include "big_float_impl.hpp"
    namespace Arageli
    {
        typedef Arageli::big_float_t<Arageli::big_float_impl> big_float;
    }
#endif

