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

template <typename P1, typename P2>
bool test_resultant_sylvester_1 (const P1& pol1, const P2& pol2)
{
    bool is_ok = true;

    typedef typename P1::coef_type Coef_type;   // WARNING! CHOOSE BETWEEN P1::coef_type and P2::coef_type.

    Coef_type resA = resultant(pol1, pol2);
    Coef_type resB = resultant(pol2, pol1);
    Coef_type resC = det(sylvester(pol1, pol2));
    Coef_type resD = det(sylvester(pol2, pol1));

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

    return is_ok;
}

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
        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("2*x^2+x+1"),
                sparse_polynom<big_int>("4*x+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("4*x^2+x+1"),
                sparse_polynom<big_int>("2*x+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("2*x^2+x+1"),
                sparse_polynom<big_int>("4*x^2+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("4*x^2+x+1"),
                sparse_polynom<big_int>("4*x+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("4*x^2+x+1"),
                sparse_polynom<big_int>("4*x^2+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("2*x^2+x+2"),
                sparse_polynom<big_int>("4*x+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("4*x^2+x+2"),
                sparse_polynom<big_int>("2*x+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("2*x^2+x+2"),
                sparse_polynom<big_int>("4*x^2+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("4*x^2+x+2"),
                sparse_polynom<big_int>("4*x+1")
            );

        is_ok &=
            test_resultant_sylvester_1
            (
                sparse_polynom<big_int>("4*x^2+x+2"),
                sparse_polynom<big_int>("4*x^2+1")
            );
    }
    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_END

    return is_ok ? resOK : resFAIL;
}
