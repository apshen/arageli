/*****************************************************************************

    hexacosichoron.hpp

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

#ifndef _ARAGELI_APP_POLYHEDRON_hexacosichoron_hpp_
#define _ARAGELI_APP_POLYHEDRON_hexacosichoron_hpp_

#include "driver.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


// 600-cell, edge = 1/(golden ratio)
class Hexacosichoron : public Processor
{
    virtual void generate (std::ostream& out, CmdArgs& cmdargs) const;
};


}
}
}

#endif
