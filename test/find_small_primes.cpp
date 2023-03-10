/*****************************************************************************

    test/find_small_primes.cpp

    This file is a part of the Arageli library.

    Copyright (C) 1999--2007 Nikolai Yu. Zolotykh
    Copyright (C) 2006 Aleksey Bader

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

#include <arageli/arageli.hpp>
#include "test_common.hpp"

TEST_FUNCTION(find_small_primes, "Test find_small_primes fuction.")
{
    using namespace Arageli;
    TestResult res = resOK;

    int N = 50;
    int *int_array = new int[N];
    small_primes(int_array, N);
    int k = 0;

    for(int i = 0; i < N && k < 5; ++i)
    {
        //std::cout << int_array[i] << ' ';
        if (!is_prime(int_array[i])) k++;
    }
    if (k > 0)
    {
        res = resFAIL;
        tout << "Test fail on " << k << "tests." << '\n';
    }
    return res;
}


/* End of file find_small_primes.cpp */
