/*****************************************************************************

    hexacosichoron.cpp

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
#include "hexacosichoron.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


void Hexacosichoron::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    if(cmdargs.dim.isSet() && cmdargs.dim.getValue() != 4)
    {
        throw "This type of polyhedron is supported for dimension 4 only.";
    }
    int dim = cmdargs.dim.getValue();
    int r = cmdargs.random.getValue();
    if(r > 0)
    {
        throw "Random generation is not supported for this type";
    }
    else
    {
        unsigned int numVertices = 120;
        out << numVertices << ' ' << 5 << '\n';

        out << std::setprecision(15);

        for ( unsigned int i = 0; i < 4; ++i )
        {
            out << "1 ";
            for ( unsigned int j = 0; j < 4; ++j )
            {
                if (i == j)
                    out << "  1";
                else
                    out << "  0";
            }
            out << '\n';

            out << "1 ";
            for ( unsigned int j = 0; j < 4; ++j )
            {
                if (i == j)
                    out << " -1";
                else
                    out << "  0";
            }
            out << '\n';
        }

        for (unsigned int i = 0; i < 16; ++i)
        {
            out << "1 ";
            for (unsigned int j = 0; j < 4; ++j)
            {
                if(i & (1 << j))
                    out << "  0.5";
                else
                    out << " -0.5";
            }
            out << '\n';
        }

        double phi = (1 + std::sqrt(5.0))*0.5;
        double phiInv = 1/phi;

        /***** 0.5*(±1, ±phi, ±phiInv, 0) even permutations*****/

        double val1[2] = {0.5, -0.5};
        double val2[2] = {phi*0.5, -phi*0.5};
        double val3[2] = {phiInv*0.5, -phiInv*0.5};
        double val4[1] = {0};

        out << std::setprecision(15);

        for (unsigned int t = 0; t < 1; ++t)
        {
            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < 2; ++j)
                {
                    for (unsigned int k = 0; k < 2; ++k)
                    {
                        out << "1 " << val1[i] << ' ' << val2[j] << ' ' << val3[k] << ' ' << val4[t] << '\n';
                        out << "1 " << val2[j] << ' ' << val3[k] << ' ' << val1[i] << ' ' << val4[t] << '\n';
                        out << "1 " << val3[k] << ' ' << val1[i] << ' ' << val2[j] << ' ' << val4[t] << '\n';
                        out << "1 " << val1[i] << ' ' << val3[k] << ' ' << val4[t] << ' ' << val2[j] << '\n';
                        out << "1 " << val2[j] << ' ' << val1[i] << ' ' << val4[t] << ' ' << val3[k] << '\n';
                        out << "1 " << val3[k] << ' ' << val2[j] << ' ' << val4[t] << ' ' << val1[i] << '\n';
                        out << "1 " << val1[i] << ' ' << val4[t] << ' ' << val2[j] << ' ' << val3[k] << '\n';
                        out << "1 " << val2[j] << ' ' << val4[t] << ' ' << val3[k] << ' ' << val1[i] << '\n';
                        out << "1 " << val3[k] << ' ' << val4[t] << ' ' << val1[i] << ' ' << val2[j] << '\n';
                        out << "1 " << val4[t] << ' ' << val3[k] << ' ' << val2[j] << ' ' << val1[i] << '\n';
                        out << "1 " << val4[t] << ' ' << val1[i] << ' ' << val3[k] << ' ' << val2[j] << '\n';
                        out << "1 " << val4[t] << ' ' << val2[j] << ' ' << val1[i] << ' ' << val3[k] << '\n';
                    }
                }
            }
        }
    }
}


}
}
}
