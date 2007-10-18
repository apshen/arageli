/*****************************************************************************

    hgcd.hpp

    This file is a part of the Arageli library.

    Copyright (C) 2007 Sergey V. Lobanov

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
    \file hgcd.hpp
    \brief <!--ADD BRIEF HEADER DESCRIPTION HERE-->

    <!--ADD ADDITIONAL FILE DESCRIPTION HERE-->
*/


#ifndef _ARAGELI__hgcd_hpp_
#define _ARAGELI__hgcd_hpp_

#include "config.hpp"

// REFERENCE ADDITIONAL HEADERS HERE


namespace Arageli
{

//x and y can be changed
template<typename X,typename Y,typename M>
void emgcd_0(X& x, Y& y,M& u1,M& u2,M& v1,M& v2,M& w1,M& w2);

//x and y can be changed
template<typename X,typename Y,typename M>
void emgcd_0(X& x, Y& y,M& u1,M& u2);//truncated version

template<typename P0,typename P1,typename M>
void gcd_emgcd_0(P0& p0,P1& p1,M& p2);//,M& p3,M& v1,M& v2,M& w1,M& w2)

//it's wrapper. How to return value by reference?
//
template<typename P>
P gcd_emgcd(P p0,P p1);

// PLACE ALL DECLARATIONS AND INLINE IMPLEMENTATIONS HERE


} // namespace Arageli


#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE__hgcd
    #include "hgcd.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE__hgcd
#endif

#endif    // #ifndef _ARAGELI__pattern_hpp_
