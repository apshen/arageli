/*****************************************************************************

    test/bug_2655678_invalid_resultant.cpp

    This file is a part of the Arageli library test base.

    Copyright (C) 2010 Sergey S. Lyalin
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
    \file bug_2655678_invalid_resultant.cpp
    \brief This file includes test for bug #2655678.

    To test bug #2655678 "Invalid work of resultant() function" and related
    stuff.
*/


#include "stdafx.hpp"

using namespace Arageli;


namespace
{

// PLACE AUXILIARY CODE HERE

}


TEST_FUNCTION
(
    resultant_bug_2655678,
    "Test bug #2655678 \"Invalid work of resultant() function.\""
)
{
    bool is_ok = true;

    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_BEGIN
    {
        sparse_polynom<big_int> pol1 = "2*x^2+x+1";
        sparse_polynom<big_int> pol2 = "4*x+1";

        big_int resA = resultant(pol1, pol2);
        big_int resB = resultant(pol2, pol1);
        big_int resC = det(sylvester(pol1, pol2));
        big_int resD = det(sylvester(pol2, pol1));

        if(resC != resD)
        {
            is_ok = false;
            tout
                << "det(sylvester(" << pol1 << ", " << pol2 << ")) == " << resC << " != "
                << resD << " == det(sylvester(" << pol2 << ", " << pol1 << "))\n";
        }

        if(resA != resB)
        {
            is_ok = false;
            tout
                << "resultant(" << pol1 << ", " << pol2 << ") == " << resA << " != "
                << resB << " == resultant(" << pol2 << ", " << pol1 << ")\n";
        }

        if(resA != resC)
        {
            is_ok = false;
            tout
                << "resultant(" << pol1 << ", " << pol2 << ") == " << resA << " != "
                << resC << " == det(sylvester(" << pol1 << ", " << pol2 << "))\n";
        }

        if(resB != resC)
        {
            is_ok = false;
            tout
                << "resultant(" << pol2 << ", " << pol1 << ") == " << resB << " != "
                << resC << " == det(sylvester(" << pol1 << ", " << pol2 << "))\n";
        }
    }
    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_END

    return is_ok ? resOK : resFAIL;
}
