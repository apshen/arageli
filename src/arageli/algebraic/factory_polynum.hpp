/*****************************************************************************

    algebraic/factory_polynum.hpp

    This file is a part of the Arageli library.

    Copyright (C) 2010 Natalia Klemenova

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


#ifndef FACTORY_POLYNUM_H
#define FACTORY_POLYNUM_H

#include "polynomial_number.hpp"
#include "../factory.hpp"

//using namespace Arageli;
namespace Arageli
{

/// Specialization of common factory template for algebraic number (PolynomialNumber)
template <>
struct factory<PolynomialNumber >
{
private:

    typedef PolynomialNumber T;

public:

    static const bool is_specialized = true;

    static const T& unit ()
    {
        static const T unit_s = T(rational<big_int>(1, 1));
        return unit_s;
    }

    static const T& unit (const T& x)
    {
        //static const T unit_s(*(x.BasisPol), rational<big_int>(1, 1)); // I
        //static const T unit_s(x);
        return unit(); //unit_s; // I
    }

    static const T& opposite_unit ()
    {
        static const T opposite_unit_s = T(rational<big_int>(-1, 1));
        return opposite_unit_s;
    }

    static const T& opposite_unit (const T& x)
    {
        return opposite_unit();
    }

    static const T& null ()
    {
        static const T null_s = T(rational<big_int>(0, 1));
        return null_s;
    }


    static const T& null (const T& x)
    {
        static const T null_s(*(x.BasisPol), rational<big_int>(0, 1)); // I
        return null(); //null_s; // I
    }

};

//
//
//inline PolynomialNumber abs(PolynomialNumber & x)
//{
//    return PolynomialNumber::abs(x); //PolynomialNumber::
//}
//



/// Specialization of common factory template for basis_field in algebraic number
/*template <>
struct factory<basis_field >
{
private:

    typedef basis_field T;

public:

    static const bool is_specialized = true;

    static const T& unit ()
    {
        static const T unit_s = T(sparse_polynom<big_int>::monom(1, 0));
        return unit_s;
    }

    static const T& unit (const T& x)
    {
        return unit();
    }

    static const T& opposite_unit ()
    {
        static const T opposite_unit_s = T(sparse_polynom<big_int>::monom(-1, 0));
        return opposite_unit_s;
    }

    static const T& opposite_unit (const T& x)
    {
        return opposite_unit();
    }

    static const T& null ()
    {
        static const T null_s = T("0");
        return null_s;
    }

    static const T& null (const T& x)
    {
        return null();
    }

};
*/
}
#endif