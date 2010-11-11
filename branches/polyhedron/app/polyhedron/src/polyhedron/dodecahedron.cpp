/*****************************************************************************

    dodecahedron.cpp

    This file is a part of Polyhedron Software, a generator of various classes
    of polyhedra.

    The Polyhedron Software is a part of the Arageli library.

    Copyright (C) 2010 Anastasya A. Ryzhova

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
#include "dodecahedron.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


void Dodecahedron::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    if(cmdargs.dim.isSet() && cmdargs.dim.getValue() != 3)
    {
        throw "This type of polyhedron is supported for dimension 3 only.";
    }
    int r = cmdargs.random.getValue();
    if(r > 0)
    {
        throw "Random generation is not supported for this type";
    }
    else
    {
        unsigned int numVertices = 20;

        out << numVertices << ' ' << 4 << '\n';

        double phi = (1 + std::sqrt(5.0))*0.5;
        double phiInv = 1/phi;

        double val1[2] = {phiInv, -phiInv};
        double val2[2] = {phi, -phi};

        unsigned int n = 1 << 3;

        for (unsigned int i = 0; i < n; ++i)
        {
            out << "1 ";
            for (unsigned int j = 0; j < 3; ++j)
            {
                if(i & (1 << j))
                    out << "  1";
                else
                    out << " -1";
            }
            out << '\n';
        }

        out << std::setprecision(15);

        for (unsigned int i = 0; i < 2; ++i)
        {
            for (unsigned int j = 0; j < 2; ++j)
            {
                out << "1 " << '0'     << ' ' << val1[i] << ' ' << val2[j] << '\n';
                out << "1 " << val1[i] << ' ' << val2[j] << ' ' << '0'     << '\n';
                out << "1 " << val2[j] << ' ' << '0'     << ' ' << val1[i] << '\n';
            }
        }
    }
}


}
}
}
