/*****************************************************************************

    cube.cpp

    This file is a part of Polyhedron Software, a generator of various classes
    of polyhedra.

    The Polyhedron Software is a part of the Arageli library.

    Copyright (C) 2010 Sergey S. Lyalin, Anastasya A. Ryzhova

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
#include "cube.hpp"
#include "utility.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


void Cube01::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    int dim = cmdargs.dim.getValue();
    int r = cmdargs.random.getValue();
    if(r > 0)
    {
        matrix<float> res(r, dim+1);

        typedef set::grid1<double> Set;
        int m = 1000000;
        Set set(0, 1, 1.0/m);
        typedef rnd::equiprob<Set> Rnd;
        Rnd rnd(set);

        for(int i = 0; i < r; ++i)
        {
            res(i, 0) = 1;
            for(int j = 1; j <= dim; ++j)
                res(i, j) = rnd();
        }

        output_matrix(out, res);
    }
    else
    {
        out << 2*dim << ' ' << dim+1 << '\n';
        for(int i = 0; i < dim; ++i)
        {
            out << "0  ";
            for(int j = 0; j < dim; ++j)
            {
                if(i == j)
                    out << " 1";
                else
                    out << " 0";
                if(j != dim - 1)
                    out << ' ';
            }
            out << '\n';

            out << "1  ";
            for(int j = 0; j < dim; ++j)
            {
                if(i == j)
                    out << "-1";
                else
                    out << " 0";
                if(j != dim - 1)
                    out << ' ';
            }
            out << '\n';
        }
    }
}

void Cube::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    int dim = cmdargs.dim.getValue();
    int r = cmdargs.random.getValue();
    if(r > 0)
    {
        throw "Random generation is not supported for this type";
    }
    else
    {
        unsigned int numVertices = 1 << dim;
        out << numVertices << ' ' << dim+1 << '\n';
        for (unsigned int i = 0; i < numVertices; ++i)
        {
            out << "1 ";
            for (unsigned int j = 0; j < dim; ++j)
            {
                if(i & (1 << j))
                    out << "  1";
                else
                    out << " -1";
            }
            out << '\n';
        }
    }
}


}
}
}
