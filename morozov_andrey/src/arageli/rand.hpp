/*****************************************************************************

    rand.hpp

    This file is a part of the Arageli library.

    Copyright (C) 1999--2007 Nikolai Yu. Zolotykh
    Copyright (C) 2005 Andrey Somsikov
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
    \file rand.hpp
    Defines T rand(T maxVal) function for each type.
*/

#ifndef __ARAGELI_rand_hpp__
#define __ARAGELI_rand_hpp__

#include <cstdlib>

#include "config.hpp"
#include "exception.hpp"
#include "factory.hpp"
#include "cmp.hpp"
#include "big_int.hpp"
#include "rational.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "sparse_polynom.hpp"

namespace Arageli
{

// WARNING! The following definition is not portable
// (but correct for Win32 platform; see commented directives).
//TODO: Make it portable.
//#if defined(WIN32)
    #define ARAGELI_RAND_GENERAL ( std::rand() )
    #define ARAGELI_RAND_SHORT ( (unsigned short)ARAGELI_RAND_GENERAL )
    #define ARAGELI_RAND_INT ( ((unsigned int)ARAGELI_RAND_SHORT << 16) | ((unsigned int)ARAGELI_RAND_SHORT) )
    #define ARAGELI_RAND_LONG ( (unsigned long)ARAGELI_RAND_INT )
    #define ARAGELI_RAND_LONG_LONG ( (((unsigned long long)ARAGELI_RAND_LONG) << 32) | ((unsigned long long)ARAGELI_RAND_LONG) )
//#else
//    #error "ARAGELI_RAND_XXX macros not defined for this platform"
//#endif

#define ARAGELI_DEFINE_FUNCTION_T_RAND_T(type, rand_value)    \
    inline type rand(type maxVal) { return 0==maxVal ? 0 : (((type)rand_value)%(maxVal+1)); }

ARAGELI_DEFINE_FUNCTION_T_RAND_T(char, ARAGELI_RAND_GENERAL);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(unsigned char, ARAGELI_RAND_GENERAL);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(short, ARAGELI_RAND_SHORT);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(unsigned short, ARAGELI_RAND_SHORT);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(int, ARAGELI_RAND_INT);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(unsigned int, ARAGELI_RAND_INT);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(long, ARAGELI_RAND_LONG);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(unsigned long, ARAGELI_RAND_LONG);
#if defined(ARAGELI_LONG_LONG_SUPPORT)
ARAGELI_DEFINE_FUNCTION_T_RAND_T(long long, ARAGELI_RAND_LONG_LONG);
ARAGELI_DEFINE_FUNCTION_T_RAND_T(unsigned long long, ARAGELI_RAND_LONG_LONG);
#endif // ARAGELI_LONG_LONG_SUPPORT

inline big_int rand(big_int maxVal)
{
    return big_int::random_in_range(maxVal);
}

template <typename T>
inline rational<T> rand(const rational<T> &maxVal)
{
    if (is_unit(maxVal.denominator()))
        return rational<T>(rand(maxVal.numerator()));

    T q(rand(maxVal.denominator()));
    if (is_null(q))
        return factory<rational<T> >::null();
    return rational<T>(rand(maxVal.numerator()), q);
}

template <typename F, typename I>
inline monom<F, I> rand(const monom<F, I> &maxVal)
{
    return monom<F, I>(rand(maxVal.coef()), rand(maxVal.degree()));
}

template <typename F, typename I, bool REFCNT>
inline sparse_polynom<F, I, REFCNT>
rand(const sparse_polynom<F, I, REFCNT> &maxVal)
{
    sparse_polynom<F, I, REFCNT> res;
    for (typename sparse_polynom<F, I, REFCNT>::monom_const_iterator it = maxVal.monoms_begin();
        it != maxVal.monoms_end(); it++)
        res += rand(*it);
    return res;
}

template <typename T, bool REFCNT>
inline matrix<T, REFCNT> rand(const matrix<T, REFCNT> &maxVal)
{
    matrix<T, REFCNT> res(maxVal.rows(), maxVal.cols());
    for (typename matrix<T, REFCNT>::size_type i = 0; i < res.rows(); i++)
        for (typename matrix<T, REFCNT>::size_type j = 0; j < res.cols(); j++)
            res[i] = rand(maxVal.el(i, j));
    return res;
}

} // namespace Arageli

#if 0
#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_rand
    #include "rand.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_rand
#endif
#endif

#endif
