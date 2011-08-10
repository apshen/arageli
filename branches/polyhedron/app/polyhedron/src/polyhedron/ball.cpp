/*****************************************************************************

    ball.cpp

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
#include "ball.hpp"
#include "utility.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


void Ball::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    int dim = cmdargs.dim.getValue();
    int r = cmdargs.random.getValue();
    if(r > 0)
    {
        typedef set::grid1<double> Set;
        int m = 1000000;
        Set set(-1, 1, 1.0/m);
        typedef rnd::equiprob<Set> Rnd;
        Rnd rnd(set);

        out << r << ' ' << dim+1 << '\n';

        out << std::setprecision(6);

        vector<double> val(dim);
        unsigned int i = 0;
        while (i < r)
        {
            for (unsigned int j = 0; j < dim; ++j)
            {
                val[j] = rnd();
            }

            if (dotsquare(val) < 1.0)
            {
                out << '1';
                for (unsigned int j = 0; j < dim; ++j)
                {
                    out << ' ' << val[j];
                }
                out << '\n';

                ++i;
            }
        }
    }
    else
    {
        throw "Error! Random key is absent in command line! Type Ball needs defined random key.";
    }
}


}
}
}