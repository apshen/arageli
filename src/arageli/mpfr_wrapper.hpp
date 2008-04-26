/*****************************************************************************
 *
 * mpfr_wrapper.hpp
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
 *     \file mpfr_wrapper.hpp
 *     \brief mpfr_t wrapper implementation
 *
 *     This module implements a C++ wrapper for MPFR object mpfr_t
 *
 **/

#ifndef _ARAGELI_mpfr_wrapper_h_
#define _ARAGELI_mpfr_wrapper_h_

#include <algorithm> //for std::max
#include <string> //for basic_string
#include <utility> //for pair
#include <mpfr.h> 
//#include "config.hpp"

namespace Arageli
{

class mpfr_wrapper
{
public:

    typedef mp_prec_t prec_t;
    typedef mp_exp_t exp_t;
    typedef mpfr_t rep_t;

    typedef enum
    {
        //IEEE Standard for Radix-Independent Floating-Point Arithmetic
        round_to_nearest = GMP_RNDN,
        round_toward_p_inf = GMP_RNDU,
        round_toward_m_inf = GMP_RNDD,
        round_toward_zero = GMP_RNDZ
    } rounding_mode_t; 

    static const prec_t max_prec = MPFR_PREC_MAX;
    static const prec_t min_prec = MPFR_PREC_MIN;

    //constructors
    mpfr_wrapper () :
        mode (get_default_rounding_mode())
    {
        std::cout << "Hello World\n";
        mpfr_init (rep);
    }

    mpfr_wrapper (prec_t prec, rounding_mode_t m) : 
        mode (m)
    {
        mpfr_init2 (rep, prec);
    }

    mpfr_wrapper (prec_t prec, unsigned long s, exp_t e, rounding_mode_t mode)
    {
        mpfr_init2(rep, prec);
        mpfr_set_ui_2exp(rep, s, e, static_cast<mp_rnd_t>(mode));
    }
 
    mpfr_wrapper (prec_t prec, long s, exp_t e, rounding_mode_t mode)
    {
        mpfr_init2(rep, prec);
        mpfr_set_si_2exp(rep, s, e, static_cast<mp_rnd_t>(mode));
    }
 
    mpfr_wrapper (char x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_si(rep, static_cast<long>(x), static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (unsigned char x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_ui(rep, static_cast<unsigned long>(x), static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (short x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_si(rep, static_cast<long>(x), static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (unsigned short x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_ui(rep, static_cast<unsigned long>(x), static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (int x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_si(rep, static_cast<long>(x), static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (unsigned int x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_ui(rep, static_cast<unsigned long>(x), static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (long x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_si(rep, x, static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (unsigned long x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_ui(rep, x, static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (const char *str, prec_t prec = get_default_precision(), rounding_mode_t m = get_default_rounding_mode(), int base = 10) :
         mode (m)
    {
        mpfr_init2(rep, prec);
        mpfr_set_str (rep, str, base, static_cast<mpfr_rnd_t>(m));
    }

    mpfr_wrapper (float x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_d(rep, static_cast<double>(x), static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (double x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_d(rep, x, static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (long double x) :
        mode (get_default_rounding_mode())
    {
        mpfr_init_set_ld(rep, x, static_cast<mpfr_rnd_t>(mode));
    }

    mpfr_wrapper (const mpfr_wrapper &b) :
        mode (b.get_rounding_mode())
    {
        mpfr_init2(rep, mpfr_get_prec(b.rep));
        mpfr_set (rep, b.rep, static_cast<mpfr_rnd_t>(b.get_rounding_mode()));
    }

    ~mpfr_wrapper ()
    {
        mpfr_clear(rep);
    }
    
    //precision and round mode
    static prec_t get_default_precision (void)
    {
        return mpfr_get_default_prec();
    }

    static void set_default_precision (long p)
    {
        return mpfr_set_default_prec (p);
    }

    static rounding_mode_t get_default_rounding_mode (void)
    {
       return static_cast<rounding_mode_t> (mpfr_get_default_rounding_mode());
    }

    static void set_default_rounding_mode ( rounding_mode_t m ) 
    {
        return mpfr_set_default_rounding_mode (static_cast<mpfr_rnd_t>(m));
    }

    void set_precision(prec_t prec)
    {
        mpfr_prec_round(rep, prec, static_cast<mp_rnd_t>(mode));
    }

    prec_t get_precision() const
    {
        return mpfr_get_prec (rep);
    }

    rounding_mode_t get_rounding_mode (void) const
    {
        return mode;
    }

    void set_rounding_mode ( rounding_mode_t m = mpfr_wrapper::round_to_nearest )
    {
         mode = m;
    }

    const rep_t & get_rep (void) const
    {
        return rep;
    }

    rep_t & get_rep (void)
    {
        return rep;
    }

    exp_t get_exp (void) const
    {
         return mpfr_get_exp(rep);
    }

    int sign() const
    {
        return mpfr_sgn(rep);
    }

    //exception related functions
    static exp_t get_emin()
    {
        return mpfr_get_emin();
    }

    static exp_t get_emax()
    {
        return mpfr_get_emax();
    }

    static bool set_emin(exp_t exp)
    {
        return mpfr_set_emin(exp);
    }

    static bool set_emax(exp_t exp)
    {
        return mpfr_set_emax(exp);
    }

    exp_t get_emin_min()
    {
        return mpfr_get_emin_min();
    }

    exp_t get_emin_max()
    {
        return mpfr_get_emin_max();
    }

    exp_t get_emax_min()
    {
        return mpfr_get_emax_min();
    }

    exp_t get_emax_max()
    {
        return mpfr_get_emax_max();
    }

    //operator =
    mpfr_wrapper & operator = (const mpfr_wrapper &other)
    {
      if (this != &other)
      {
           set_rounding_mode(other.get_rounding_mode());
           set_precision(other.get_precision());
           mpfr_set(rep, other.rep, static_cast<mpfr_rnd_t>(mode));
      }
      return *this;
    }

    mpfr_wrapper & operator = (char x)
    {
        mpfr_set_si(rep, static_cast<long>(x), static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (unsigned char x)
    {
        mpfr_set_ui(rep, static_cast<unsigned long>(x), static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (short x)
    {
        mpfr_set_si(rep, static_cast<long>(x), static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (unsigned short x)
    {
        mpfr_set_ui(rep, static_cast<unsigned long>(x), static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (int x)
    {
        mpfr_set_si(rep, static_cast<long>(x), static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (unsigned int x)
    {
        mpfr_set_ui(rep, static_cast<unsigned long>(x), static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (long x)
    {
        mpfr_set_si(rep, x, static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (unsigned long x)
    {
        mpfr_set_ui(rep, x, static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (float x)
    {
        mpfr_set_d(rep, static_cast<double>(x), static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (double x)
    {
        mpfr_set_d(rep, x, static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (long double x)
    {
        mpfr_set_ld(rep, x, static_cast<mpfr_rnd_t>(mode));
        return *this;
    }

    mpfr_wrapper & operator = (const char *str)
    {
        mpfr_init_set_str (rep, str, 10, static_cast<mpfr_rnd_t>(mode));
    }
    //conversion to native types 
    operator char () const
    {
        return static_cast<char> (mpfr_get_si(get_rep(), get_mp_rounding_mode()));
    }

    operator short () const
    {
        return static_cast<short> (mpfr_get_si(get_rep(), get_mp_rounding_mode()));
    }

    operator int () const
    {
        return static_cast<int> (mpfr_get_si(get_rep(), get_mp_rounding_mode()));
    }

    operator long () const
    {
        return mpfr_get_si(get_rep(), get_mp_rounding_mode());
    }

    operator unsigned char () const
    {
        return static_cast<unsigned char> (mpfr_get_ui(get_rep(), get_mp_rounding_mode()));
    }

    operator unsigned short () const
    {
        return static_cast<unsigned short> (mpfr_get_ui(get_rep(), get_mp_rounding_mode()));
    }

    operator unsigned int () const
    {
        return static_cast<unsigned int> (mpfr_get_si(get_rep(), get_mp_rounding_mode()));
    }

    operator unsigned long () const
    {
        return mpfr_get_ui(get_rep(), get_mp_rounding_mode());
    }

    operator float () const
    {
        return static_cast<float>(mpfr_get_d(get_rep(), get_mp_rounding_mode()));
    }

    operator double () const
    {
        return mpfr_get_d(get_rep(), get_mp_rounding_mode());
    }

    operator long double () const
    {
        return mpfr_get_ld(get_rep(), get_mp_rounding_mode());
    }

    mp_rnd_t get_mp_rounding_mode() const
    {
        return static_cast<mp_rnd_t>(mode);
    }

    bool issingular () const
    {
        return mpfr_nan_p (get_rep()) || mpfr_inf_p(get_rep()) || mpfr_zero_p (get_rep());
    }

    //special numbers
    bool isnan() const
    {
        return mpfr_nan_p(get_rep());
    }

    bool isinf() const
    {
        return mpfr_inf_p(get_rep());
    }

    bool isminf() const
    {
        return isinf() && static_cast<bool>(mpfr_signbit(get_rep()));
    }

    bool ispinf() const
    {
        return isinf() && !static_cast<bool>(mpfr_signbit(get_rep()));
    }

    bool isfinite() const
    {
        return mpfr_number_p(get_rep());
    }

    bool iszero() const
    {
        return mpfr_zero_p(get_rep()); 
    }

    bool ismzero() const
    {
        return iszero() && static_cast<bool>(mpfr_signbit(get_rep()));
    }

    bool ispzero() const
    {
        return iszero() && !static_cast<bool>(mpfr_signbit(get_rep()));
    }

    std::pair<std::basic_string<char>, exp_t> to_string (std::size_t n, int base) const
    {
        exp_t exp;
        char *cm = mpfr_get_str (0, &exp, base, n, get_rep(), get_mp_rounding_mode());
        std::basic_string<char> man(cm);
        mpfr_free_str(cm);
        return std::make_pair(man, exp);
    }

private:

    mpfr_t rep;
    rounding_mode_t mode;

/**/
};

#define _ARAGELI_MPFR_WRAPPER_STANDART_OPERATION(OP) \
inline mpfr_wrapper OP (const mpfr_wrapper &b, const mpfr_wrapper &c, mpfr_wrapper::prec_t prec, mpfr_wrapper::rounding_mode_t m)\
{\
    mpfr_wrapper ret(prec, m);\
    mpfr_##OP (ret.get_rep(), b.get_rep(), c.get_rep(), static_cast<mpfr_rnd_t>(m));\
    return ret;\
}

_ARAGELI_MPFR_WRAPPER_STANDART_OPERATION(add)
_ARAGELI_MPFR_WRAPPER_STANDART_OPERATION(sub)
_ARAGELI_MPFR_WRAPPER_STANDART_OPERATION(mul)
_ARAGELI_MPFR_WRAPPER_STANDART_OPERATION(div)

#undef _ARAGELI_MPFR_WRAPPER_STANDART_OPEARTION

//special functions

#define _ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(FUNC)\
inline mpfr_wrapper FUNC (const mpfr_wrapper &b, mpfr_wrapper::prec_t prec, mpfr_wrapper::rounding_mode_t m)\
{\
    mpfr_wrapper ret(prec, m);\
    mpfr_##FUNC (ret.get_rep(), b.get_rep(), static_cast<mp_rnd_t>(m));\
    return ret;\
}

#define _ARAGELI_MPFR_WRAPPER_SPECIAL_BFUNC(FUNC)\
inline mpfr_wrapper FUNC (const mpfr_wrapper &b, const mpfr_wrapper &c, mpfr_wrapper::prec_t prec, mpfr_wrapper::rounding_mode_t m)\
{\
    mpfr_wrapper ret (prec,m);\
    mpfr_##FUNC(ret.get_rep(),b.get_rep(), c.get_rep(), static_cast<mp_rnd_t>(m));\
    return ret;\
}

#define _ARAGELI_MPFR_WRAPPER_SPECIAL_MIXED_BFUNC(FUNC)\
inline mpfr_wrapper FUNC (long x, const mpfr_wrapper &b, mpfr_wrapper::prec_t prec, mpfr_wrapper::rounding_mode_t m)\
{\
    mpfr_wrapper ret (prec,m);\
    mpfr_##FUNC(ret.get_rep(), x, b.get_rep(), static_cast<mp_rnd_t>(m));\
    return ret;\
}

#define _ARAGELI_MPFR_WRAPPER_SPECIAL_VFUNC(FUNC)\
inline mpfr_wrapper FUNC (mpfr_wrapper::prec_t prec, mpfr_wrapper::rounding_mode_t m)\
{\
    mpfr_wrapper ret (prec,m);\
    mpfr_const_##FUNC(ret.get_rep(),static_cast<mp_rnd_t>(m));\
    return ret;\
}

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(log)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(log2)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(log10)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(exp)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(exp2)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(exp10)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(sin)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(cos)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(tan)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(sec)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(csc)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(cot)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(asin)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(acos)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(atan)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(sinh)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(cosh)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(tanh)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(sech)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(csch)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(coth)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(asinh)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(acosh)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(atanh)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(log1p)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(expm1)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(eint)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(erf)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(erfc)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(gamma)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(lngamma)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(zeta)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(j0)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(j1)
_ARAGELI_MPFR_WRAPPER_SPECIAL_MIXED_BFUNC(jn)

_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(y0)
_ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC(y1)
_ARAGELI_MPFR_WRAPPER_SPECIAL_MIXED_BFUNC(yn)

_ARAGELI_MPFR_WRAPPER_SPECIAL_VFUNC(pi)
_ARAGELI_MPFR_WRAPPER_SPECIAL_VFUNC(log2)
_ARAGELI_MPFR_WRAPPER_SPECIAL_VFUNC(euler)
_ARAGELI_MPFR_WRAPPER_SPECIAL_VFUNC(catalan)

_ARAGELI_MPFR_WRAPPER_SPECIAL_BFUNC(atan2)
_ARAGELI_MPFR_WRAPPER_SPECIAL_BFUNC(hypot)

#undef _ARAGELI_MPFR_WRAPPER_SPECIAL_UFUNC
#undef _ARAGELI_MPFR_WRAPPER_SPECIAL_BFUNC
#undef _ARAGELI_MPFR_WRAPPER_SPECIAL_VFUNC
#undef _ARAGELI_MPFR_WRAPPER_SPECIAL_MIXED_BFUNC

inline mpfr_wrapper fac (unsigned long x, mpfr_wrapper::prec_t prec, mpfr_wrapper::rounding_mode_t m)
{
    mpfr_wrapper ret(prec, m);
    mpfr_fac_ui(ret.get_rep(), x, static_cast<mp_rnd_t>(m));
    return ret;
}

inline mpfr_wrapper zeta (unsigned long x, mpfr_wrapper::prec_t prec, mpfr_wrapper::rounding_mode_t m)
{
    mpfr_wrapper ret(prec, m);
    mpfr_zeta_ui(ret.get_rep(), x, static_cast<mp_rnd_t>(m));
    return ret;
}

//comparison TODO add NaN handling! RESOLVED NaN will be handled in comparison operators (e.g ==. !=, >, <, etc)
inline int cmp(const mpfr_wrapper &b, const mpfr_wrapper &c)
{
    return mpfr_cmp(b.get_rep(), c.get_rep());
}

#define _ARAGELI_MPFR_WRAPPER_SIGNED_COMPARISON(NATIVE)\
inline int cmp(const mpfr_wrapper &b, NATIVE c)\
{\
    return mpfr_cmp_si (b.get_rep(), static_cast<long>(c));\
}

#define _ARAGELI_MPFR_WRAPPER_UNSIGNED_COMPARISON(NATIVE)\
inline int cmp(const mpfr_wrapper &b, NATIVE c)\
{\
    return mpfr_cmp_ui (b.get_rep(), static_cast<unsigned long>(c));\
}

_ARAGELI_MPFR_WRAPPER_SIGNED_COMPARISON(signed char)
_ARAGELI_MPFR_WRAPPER_SIGNED_COMPARISON(short)
_ARAGELI_MPFR_WRAPPER_SIGNED_COMPARISON(int)
_ARAGELI_MPFR_WRAPPER_SIGNED_COMPARISON(long)

_ARAGELI_MPFR_WRAPPER_UNSIGNED_COMPARISON(unsigned char)
_ARAGELI_MPFR_WRAPPER_UNSIGNED_COMPARISON(unsigned short)
_ARAGELI_MPFR_WRAPPER_UNSIGNED_COMPARISON(unsigned int)
_ARAGELI_MPFR_WRAPPER_UNSIGNED_COMPARISON(unsigned long)

#undef _ARAGELI_MPFR_WRAPPER_SIGNED_COMPARISON
#undef _ARAGELI_MPFR_WRAPPER_UNSIGNED_COMPARISON

inline int cmp(const mpfr_wrapper &b, float c)
{
    return mpfr_cmp_d(b.get_rep(), static_cast<double>(c));
}

inline int cmp(const mpfr_wrapper &b, double c)
{
    return mpfr_cmp_d(b.get_rep(), c);
}

inline int cmp(const mpfr_wrapper &b, long double c)
{
    return mpfr_cmp_ld(b.get_rep(), c);
}

//input and putput
std::istream & operator >> (std::istream &os, mpfr_wrapper &b);

//integer related functions
#define _ARAGELI_MPFR_WRAPPER_INT_RELATED_FUNC(FUNC)\
inline mpfr_wrapper FUNC (const mpfr_wrapper &b)\
{\
    mpfr_wrapper ret (b.get_precision(), b.get_rounding_mode());\
    mpfr_##FUNC(ret.get_rep(), b.get_rep());\
    return ret;\
}

_ARAGELI_MPFR_WRAPPER_INT_RELATED_FUNC(ceil)
_ARAGELI_MPFR_WRAPPER_INT_RELATED_FUNC(floor)
_ARAGELI_MPFR_WRAPPER_INT_RELATED_FUNC(round)
_ARAGELI_MPFR_WRAPPER_INT_RELATED_FUNC(trunc)

#undef _ARAGELI_MPFR_WRAPPER_INT_RELATED_FUNC

inline mpfr_wrapper rint (const mpfr_wrapper &b)
{
    mpfr_wrapper ret (b.get_precision(), b.get_rounding_mode());
    mpfr_rint(ret.get_rep(), b.get_rep(), b.get_mp_rounding_mode());
    return ret;
}

inline mpfr_wrapper frac (const mpfr_wrapper &b)
{
    mpfr_wrapper ret (b.get_precision(), b.get_rounding_mode());
    mpfr_frac(ret.get_rep(), b.get_rep(), b.get_mp_rounding_mode());
    return ret;
}

//miscellaneous functions
inline mpfr_wrapper nexttoward (const mpfr_wrapper &b, const mpfr_wrapper &c)
{
    mpfr_wrapper ret(b);
    mpfr_nexttoward (ret.get_rep(), c.get_rep());
    return ret;
}

inline mpfr_wrapper nextabove(const mpfr_wrapper &b)
{
    mpfr_wrapper ret(b);
    mpfr_nextabove(ret.get_rep());
    return ret;
}

inline mpfr_wrapper nextbelow(const mpfr_wrapper &b)
{
    mpfr_wrapper ret(b);
    mpfr_nextbelow(ret.get_rep());
    return ret;
}

mpfr_wrapper min (const mpfr_wrapper &b, const mpfr_wrapper &c);
mpfr_wrapper max (const mpfr_wrapper &b, const mpfr_wrapper &c);

inline int sign(const mpfr_wrapper &b)
{
    return b.sign();
}

}//namespace Arageli

//#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
//    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE__mpfr_wrapper
//    #include "mpfr_wrapper.cpp"
//    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE__mpfr_wrapper
//#endif

#endif
