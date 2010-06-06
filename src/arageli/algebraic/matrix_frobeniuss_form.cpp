/*****************************************************************************

    algebraic/matrix_frobeniuss_form.cpp

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


#include "matrix_frobeniuss_form.hpp"

namespace Arageli
{

void Frobeniuss(int degr, sparse_polynom< big_int> g, matrix<big_int> &F)
{
    int j = degr-1;
    int i = degr-1;
    int dec = 0;
    monoms mall = g.monoms_begin();

    if (mall->degree() != 0) dec = mall->degree();
    if (dec >= 1) dec++;

    for(monoms mi = g.monoms_begin(), mj = g.monoms_end(); mi != mj && mi->degree() < degr; ++mi)
    {
        while (dec > 1 && i > 0)
        {
            i--;
            F(j,i) = 1;
            j--;
            dec--;
        }

        F(0,i) = -mi->coef();
        if (i > 0) F(j,i-1) = 1;
        i--;
        j--;
        mall++;

        dec = mall->degree() - mi->degree();

        while (dec > 1 && i > 0)
        {
            i--;
            F(j,i) = 1;
            j--;
            dec--;
        }
    }
}

void Resz(int m, matrix<big_int> &F)
{
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            F(i,j) = 0;
}

} //- end namespace Arageli --------------------------------------------------------