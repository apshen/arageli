/*****************************************************************************

    utility.hpp

    This file is a part of Polyhedron Software, a generator of various classes
    of polyhedra.

    The Polyhedron Software is a part of the Arageli library.

    Copyright (C) 2010 Sergey S. Lyalin

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

#ifndef _ARAGELI_APP_POLYHEDRON_utility_hpp_
#define _ARAGELI_APP_POLYHEDRON_utility_hpp_

#if 0

namespace Arageli
{
namespace app
{
namespace polyhedron
{

typedef matrix<big_int> Matrix;

template <typename M>
M make_cone (const M& x);

template <typename M>
void output_matrix (std::ostream& out, const M& x);

Matrix cyclic_verts (int d, int n);

template <typename M>
M simplex_verts (int d);

template <typename M>
M make_verts_product (const M& x, const M& y);


}
}
}

#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_app_polyhedron_utility
    #include "utility.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_app_polyhedron_utility
#endif

#else

namespace Arageli
{
namespace app
{
namespace polyhedron
{

typedef matrix<big_int> Matrix;

Matrix cyclic_verts (int d, int n);

template <typename M>
M make_cone (const M& x)
{
    M res = x;
    res.insert_col(0, unit<typename M::value_type>());
    return res;
}

template <typename M>
void output_matrix (std::ostream& out, const M& x)
{
    out << x.nrows() << ' ' << x.ncols() << '\n';
    output_aligned(out, x, " ", " ", " ");
}

template <typename M>
M simplex_verts (int d)
{
    M res(d+1, d);

    for(int i = 0; i < d; ++i)
        res(i+1, i) = 1;

    return res;
}

template <typename M>
M make_verts_product (const M& x, const M& y)
{
    if(x.is_empty())return y;
    if(y.is_empty())return x;

    M res(x.nrows()*y.nrows(), x.ncols() + y.ncols());

    for(int i = 0; i < x.nrows(); ++i)
        for(int j = 0; j < y.nrows(); ++j)
        {
            int ij = i*y.nrows() + j;
            for(int ki = 0; ki < x.ncols(); ++ki)
                res(ij, ki) = x(i, ki);
            for(int kj = 0; kj < y.ncols(); ++kj)
                res(ij, kj + x.ncols()) = y(j, kj);
        }

    return res;
}


}
}
}

#endif

#endif