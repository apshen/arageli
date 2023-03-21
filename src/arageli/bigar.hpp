/*****************************************************************************

    bigar.hpp

    This file is a part of the Arageli library.

    WARNING. This is internal header of Arageli library.
    Do not use it directly in you code other then Arageli code.

    Copyright (C) 1999 -- 2005, 2008 Nikolai Yu. Zolotykh

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
    \file bigar.hpp
    \brief Dealing with chains of big digits as a big integer numbers.

    A lot of internal constants. Type for the big digit is being choosen.
    The four basic operations on big digits sequences (big integers).

    The implementation of the four basic arithmetic operations based on
    an assumption that there are doubledigit type and extendeddigit type
    that have double size of the digit and size greater by 1 bit then
    the digit correspondingly.
*/

#ifndef __ARAGELI_bigar_hpp__
#define __ARAGELI_bigar_hpp__

#include "config.hpp"

#include <cstddef>
#include <cstdlib>
#include <limits>

#include "std_import.hpp"


namespace Arageli { namespace _Internal
{

    // WARNING for the whole file!
    // Sizes, types and masks can faile on some platforms.
    // This is a temporary implementation.

template<
    unsigned I,
    typename std::enable_if<
        std::numeric_limits<unsigned long long int>::digits == I,
        bool>::type = true
>
unsigned long long int digit_type_detector_helper();

template<
    unsigned I,
    typename std::enable_if<
        std::numeric_limits<unsigned long int>::digits == I &&
        std::numeric_limits<unsigned long int>::digits < std::numeric_limits<unsigned long long int>::digits,
        bool>::type = true
>
unsigned long int digit_type_detector_helper();

template<
    unsigned I,
    typename std::enable_if<
        std::numeric_limits<unsigned int>::digits == I &&
        std::numeric_limits<unsigned int>::digits < std::numeric_limits<unsigned long int>::digits,
        bool>::type = true
>
unsigned int digit_type_detector_helper();

template<
    unsigned I,
    typename std::enable_if<
        std::numeric_limits<unsigned short int>::digits == I &&
        std::numeric_limits<unsigned short int>::digits < std::numeric_limits<unsigned int>::digits,
        bool>::type = true
>
unsigned short int digit_type_detector_helper();

#ifdef ARAGELI_HAS_UINT128_SUPPORT
typedef decltype(digit_type_detector_helper<64>()) digit;
typedef __uint128_t doubledigit;
#else
typedef decltype(digit_type_detector_helper<32>()) digit;
typedef decltype(digit_type_detector_helper<64>()) doubledigit;
#endif

typedef doubledigit extendeddigit;
typedef unsigned short int bit;
constexpr digit max_digit = digit(-1);
constexpr int bits_per_digit = std::numeric_limits<digit>::digits;
constexpr extendeddigit BASE = extendeddigit(1) << bits_per_digit;

std::size_t do_big_int_to_bdn
(
    digit* a,
    digit* b,
    std::size_t n,
    digit bdn_radix
);

std::size_t do_bdn_to_big_int
(
    digit* a,
    digit* b,
    std::size_t n,
    digit bdn_radix
);

digit do_add_1
(
    digit* p1,
    digit a,
    std::size_t m
);

digit do_add
(
    digit* p1,
    const digit* p2,
    std::size_t m,
    std::size_t n
);

std::size_t do_add_and_set_carry
(
    digit* p1,
    const digit* p2,
    std::size_t m,
    std::size_t n
);

int do_sub
(
    digit* p1,
    const digit* p2,
    std::size_t m,
    std::size_t n
);

std::size_t do_optimize
(
    const digit* a,
    std::size_t n
);

std::size_t do_mult_classic
(
    const digit* u,
    const digit* v,
    digit* w,
    std::size_t m,
    std::size_t n
);

std::size_t do_mult
(
    const digit* u,
    const digit* v,
    digit* w,
    std::size_t m,
    std::size_t n
);

digit do_divide_by_digit
(
    const digit* a,
    digit* p,
    std::size_t n,
    digit d
);

std::size_t do_divide
(
    digit* u,
    digit* v,
    digit* q,
    std::size_t m,
    std::size_t n
);

} // namesace _Internal
} // namespece Arageli

#endif
