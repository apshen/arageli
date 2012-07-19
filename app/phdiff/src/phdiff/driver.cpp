/*****************************************************************************

    driver.cpp

    This file is a part of phdiff tool, comparator for polyhedrons.
    Phdiff tool is a part of the Arageli library.

    Copyright (C) 2012 Anastasya A. Ryzhova

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
#include "driver.hpp"

namespace Arageli
{
namespace app
{
namespace phdiff
{


CmdArgs::CmdArgs (TCLAP::CmdLine& cmd) :
    type
    (
        "t",
        "type",

        "Type of polyhedra to compare.\n"
        "int.\n "
        "bigint.\n"
		"double.\n",
        
        true,
        "",
        "int, bigint, double",
        cmd
    ),
    epsilon
    (
        "e",
        "epsilon",

        "Accurace.\n",

        false,
        0,
        "positive float number",
        cmd
    ),
    epsilon_adjust
    (
        "a",
        "epsilon_adjust",

        "Turns automatic adjustment of epsilon on.",

        cmd
    )
	{
    }

}
}
}
