/*****************************************************************************

    sphere.cpp

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
#include "sphere.hpp"
#include "utility.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


#define HOLD_AS_MATRIX 0

void Sphere::generate (std::ostream& out, CmdArgs& cmdargs) const
{
    int dim = cmdargs.dim.getValue();
    int r = cmdargs.random.getValue();
    if(r > 0)
    {
        typedef set::grid1<double> Set;
        int m = 1000000;
        Set set(-1, 1, 1.0/m);
        typedef rnd::equiprob<Set> Rnd;
        Rnd rnd(10, set);

        out << r << ' ' << dim+1 << '\n';

        out << std::setprecision(6);

#if HOLD_AS_MATRIX

        matrix<rational<> > mat(r, dim);

#endif

        vector<double> val(dim);
        unsigned int i = 0;
        while (i < r)
        {
            for (unsigned int j = 0; j < dim; ++j)
            {
                val[j] = rnd();
            }

            double len = std::sqrt(dotsquare(val));

            if (len < 1.0)
            {
                val /= len;
                out << '1';
                for (unsigned int j = 0; j < dim; ++j)
                {
#if HOLD_AS_MATRIX

                    mat(i, j) = val[j];
#endif
                    out << ' ' << val[j];
                }
                out << '\n';

                ++i;
            }
        }

#if HOLD_AS_MATRIX

        output_aligned(std::cout, mat);
        polyhedron<> p(mat, fromvert);
        std::ofstream file("polyhedron.wrl");
        output_vrml(file, p);

#endif

    }
    else
    {
        throw "Error! Random key is absent in command line! Type Ball needs defined random key.";
    }
}


}
}
}
