/*****************************************************************************

    algebraic/matrix_frobeniuss.hpp

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


#ifndef ALGEBRAIC_matrix_frobeniuss_form_hpp
#define ALGEBRAIC_matrix_frobeniuss_form_hpp

#include "../config.hpp"

#include "../matrix.hpp"
#include "../sparse_polynom.hpp"
#include "../rational.hpp"
#include "../big_int.hpp"

namespace Arageli
{

typedef sparse_polynom< big_int>::monom_iterator monoms;
typedef sparse_polynom< rational< big_int>>::monom_iterator monoms_r;


void Frobeniuss(int degr, sparse_polynom< big_int> g, matrix<big_int> &F);

void Resz(int m, matrix<big_int> &F);

} //- end namespace Arageli --------------------------------------------------------
#endif