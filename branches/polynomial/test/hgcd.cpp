/*****************************************************************************

    test/hgcd.cpp

    This file is a part of the Arageli library test base.

    Copyright (C) 2008 Sergey V. Lobanov

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
    \file hgcd.cpp
    \brief This file includes test for hgcd.

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#include "stdafx.hpp"

using namespace Arageli;


namespace
{
template <typename P>
bool gcd_emgcd_test()
{
	//__asm int 3;
	P p0="5*x^17+7*x^16+90";
	P p1="5*x^20+x^10+1";
	P p2="2*x^18+x^6+2";
	P p3="3*x^7+x^3+3";
	p1*=p0;
	p2*=p1;
	p3*=p1;
	if (gcd_emgcd(p2,p3)!=(gcd(p2,p3)+0))
	{
		tout<<"\nFAILED:"<<typeid(P).name()<<"\np2="<<p2<<"\np3="<<p3<<"\n";
		return false;
	}

	return true;
}
// PLACE AUXILIARY CODE HERE

}




TEST_FUNCTION
(
    gcd_emgcd,
    "Test for gcd_emgcd function."
)
{
    bool is_ok = true;

    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_BEGIN
    {
        // PLACE TEST CODE HERE

        // FOR EXAMPLE:
         is_ok &= gcd_emgcd_test<sparse_polynom<rational<> > >();
         //is_ok &= gcd_emgcd_test<polynom<rational<> > >(); //compilation problem
    }
    ARAGELI_TS_ALLEXCEPT_CATCH_REGION_END

    return is_ok ? resOK : resFAIL;
}
