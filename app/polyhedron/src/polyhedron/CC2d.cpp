/*****************************************************************************

    CC2d.cpp

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

#include "stdafx.hpp"
#include "CC2d.hpp"
#include "utility.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


void CC2d::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    int dim = cmdargs.dim.getValue();
    int n = cmdargs.number.getValue();
    if(dim <= 2)
    {
        std::cerr << "[ERROR] You should provide dim >= 3.";
        return;
    }
    if(n <= 0)
    {
        std::cerr << "[ERROR] You should provide positive integer n.";
        return;
    }

    int a = intsqrt(dim);
    if(a*a < dim)++a;
    int b = dim/a;
    int c = dim%a;

    matrix<big_int> c2a = cyclic_verts(2*a, n);
    matrix<big_int> c2c = cyclic_verts(2*c, n);

    matrix<big_int> c2as;
    for(int i = 0; i < b; ++i)
    {
        c2as = make_verts_product(c2as, c2a);
    }

    c2as = make_verts_product(c2as, c2c);
    output_matrix(out, make_cone(c2as));
}


}
}
}
