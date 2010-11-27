/*****************************************************************************

    permute.hpp

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

#ifndef _ARAGELI_APP_POLYHEDRON_permute_hpp_
#define _ARAGELI_APP_POLYHEDRON_permute_hpp_

#if 0

namespace Arageli
{
namespace app
{
namespace polyhedron
{


template <typename T>
void AllPermutations(std::ostream& out, std::vector<T>& val);


}
}
}

#ifdef ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE
    #define ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_app_polyhedron_permute
    #include "permute.cpp"
    #undef  ARAGELI_INCLUDE_CPP_WITH_EXPORT_TEMPLATE_app_polyhedron_permute
#endif

#else

#include "permute.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


template <typename T>
bool IsValueRepeated(std::vector<T>& val, unsigned int ind, unsigned int j)
{
    for (unsigned int i = 0; i < j; ++i)
    {
        if (val[i+ind] == val[j+ind])
            return true;
    }

    return false;
}

template <typename T>
void Permute(std::ostream& out, std::vector<T>& val, unsigned int m, unsigned int n, unsigned int ind)
{
    std::vector<T> valTemp(val);

    if (n > 1)
    {
        Permute(out, val, m, n-1, ind+1);
        for (unsigned int i = 1; i < n; ++i)
        {
            if (!IsValueRepeated(valTemp, ind, i))
            {
                val[ind] = valTemp[i+ind];
                val[i+ind] = valTemp[ind];

                Permute(out, val, m, n-1, ind+1);

                val[ind] = valTemp[ind];
                val[i+ind] = valTemp[i+ind];
            }
        }
    }
    else
    {
        out << '1';
        for (unsigned int i = 0; i < m; ++i)
        {
            out << ' ' << val[i];
        }
        out << '\n';
    }
}

template <typename T>
void AllPermutations(std::ostream& out, std::vector<T>& val)
{
    unsigned int m = val.size();
    Permute(out, val, m, m, 0);
}


}
}
}

#endif

#endif
