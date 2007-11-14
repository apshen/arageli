/*****************************************************************************

    karatsuba_vs_classic.cpp

    This file is a part of the Arageli library.

    Copyright (C) 1999--2007 Nikolai Yu. Zolotykh
    Copyright (C) 2005--2007 Sergey S. Lyalin
    Copyright (C) 2007 Aleksey Bader
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
    \file karatsuba_vs_classic.cpp
    \brief Utility to determine threshold value.
    This value is minimum number length for which we karatsuba multiplication algorithm.
    Classic multiplication algorithm is used for numbers with length less than threshold value.

*/


#include <arageli/arageli.hpp>

using namespace Arageli;
// We need internal namespace because of using interal big_int representation.
using namespace Arageli::_Internal;

int main ()
{
    const unsigned int lower_bound = 5;
    const unsigned int upper_bound = 100;
    const unsigned int repeat = 10;
    const unsigned int tolerance = 5;
    unsigned int threshold = 0;
    timer karatsuba_timer(false), classic_timer(false);
    double karatsuba_duration = 0., classic_duration = 0.;
    do{
        big_int 
    } while ((upper_bound - lower_bound) > tolerance)
    big_int a = random_with_length(), b = random_with_length();
    return threshold;
}
