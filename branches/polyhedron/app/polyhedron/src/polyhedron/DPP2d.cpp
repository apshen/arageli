/*****************************************************************************

    DD2d.cpp

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
#include "DPP2d.hpp"
#include "utility.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


void DPP2d::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    int dim = cmdargs.dim.getValue();
    int n = cmdargs.number.getValue();
    if(dim < 2)
    {
        std::cerr << "[ERROR] You should provide dim >= 2.";
        return;
    }
    if(n < 3)
    {
        std::cerr << "[ERROR] You should provide n >= 3.";
        return;
    }

    Matrix res(2*dim + dim*(n-3) + dim + 1, 2*dim + 1);

    for(int i = 0; i < dim; ++i)
    {
        // Main inequalities.

        res(i, dim + 1 + i) = 1;
        res(i + dim, i + 1) = n;
        res(i + dim, dim + 1 + i) = -1;

        for(int k = 0; k < n-3; ++k)
        {
            res(2*dim + dim*k + i, i + 1) = -(2*k+1);
            res(2*dim + dim*k + i, dim + i + 1) = -1;
            res(2*dim + dim*k + i, 0) = (2*k+1)*(n+k) - k*k + n*n;
        }

        res(2*dim + dim*(n-3) + i, i + 1) = -(2*n-3);
        res(2*dim + dim*(n-3) + i, dim + i + 1) = -1;
        res(2*dim + dim*(n-3) + i, 0) = 2*n*(2*n-3);

        // Dwarfed inequality.

        res(res.nrows()-1, i+1) = -1;
    }

    res(res.nrows()-1, 0) = 2*n-1;

    output_matrix(out, res);
}


}
}
}
