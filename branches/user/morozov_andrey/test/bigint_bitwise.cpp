/*****************************************************************************

    test/bigint_bitwise.cpp

    This file is a part of the Arageli library test base.

    Copyright (C) 1999--2007 Nikolai Yu. Zolotykh
    Copyright (C) 2005--2007 Sergey S. Lyalin
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
    \file bigint_bitwise.cpp
    \brief This file includes test for bitwise operations with big_int.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#include "stdafx.hpp"

using namespace Arageli;


namespace
{

bool operator_and ()
{
    for(int i = 0; i < 1000; ++i)
        for(int j = 10; j < 100; ++j)
        {
            big_int
                a = big_int::random_with_length_or_less(j),
                b = big_int::random_with_length_or_less(j + i);
            big_int c = a & b;

            if(c.length() > a.length() || c.length() > b.length())
            {
                tout
                    << "FAILED 1 with:"
                    << "\n\ta = " << a
                    << "\n\tb = " << b
                    << "\n\tc = " << c
                    << "\n\ta.length() = " << a.length()
                    << "\n\tb.length() = " << b.length()
                    << "\n\tc.length() = " << c.length() << '\n';
                return false;
            }

            for(int k = 0; k < c.length(); ++k)
                if(c[k] != a[k] && b[k])
                {
                    tout
                        << "FAILED 2 with:"
                        << "\n\ta = " << a
                        << "\n\tb = " << b
                        << "\n\tc = " << c
                        << "\n\ta.length() = " << a.length()
                        << "\n\tb.length() = " << b.length()
                        << "\n\tc.length() = " << c.length()
                        << "\n\tk = " << k
                        << "\n\ta[k] = " << a[k]
                        << "\n\tb[k] = " << b[k]
                        << "\n\tc[k] = " << c[k] << '\n';
                    return false;
                }
        }

    return true;
}

}


TEST
(
    big_int,
    bitwise_operators,
    "Test for bitwise operations with big_int."
)
{
    bool is_ok = true;

    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_BEGIN
    {
        is_ok &= operator_and();
    }
    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_END

    return is_ok ? resOK : resFAIL;
}
