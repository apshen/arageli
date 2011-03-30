/*****************************************************************************

    driver.cpp

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
#include "driver.hpp"

namespace Arageli
{
namespace app
{
namespace polyhedron
{


CmdArgs::CmdArgs (TCLAP::CmdLine& cmd) :
    type
    (
        "t",
        "type",

        "Type of polyhedron to generate.\n"
        "simplex, tetrahedron -- standard d-simplex with vertices (0, 0, ..., 0), (1, 0, ..., 0), ..., (0, 0, ..., 1). "
        "Simplex is represented as a set of inequalities by the default.\n"
        "cube01 -- d-dimensional cube [0, 1]. Cube01 is represetnted as a set of inequalities by the default.\n"
        "cyclic -- d-dimensional cyclic polytop built with the cure (t, t^2, ..., t^d). Cyclic polytopes are represented by a set of vertices by the default.\n"
        "dwarfedcube, dwarfcube -- d-dimensional dwarfed cube [0, 2] with the dwarfing inequality x_1 + x_2 + ... + x_d <= 3. "
        "Dwarfed cubes are represented as a set of inequalities by the default.\n"
        "CC2d -- 2*d-dimensioanal sqrt(d)-fold product of (2*sqrt(d))-dimensional cyclic polytopes with n vertices each. "
        "This class of polyhedra is represented as a set of vertices by the default.\n"
        "DPP2d -- 2*d-dimensional dwarfed product of polygons. Polytopes from this class are represeted as a set of inequalities by the default.\n"
        "TT2d -- 2*d-dimensional product of two d-dimensional simplices. This class of polyhedra is represented as a set of vertices by the default.",

        true,
        "",
        "tetrahedron, simplex, cube01, cube, octahedron, dodecahedron, icosahedron, icositetrachoron, hecatonicosachoron, hexacosichoron, ball, sphere, cyclic, dwarfedcube, dwarfcube, CC2d, DPP2d, TT2d",
        cmd
    ),
    dim
    (
        "d",
        "dim",

        "Base dimension for polyhedron.",

        false,
        0,
        "positive integer number",
        cmd
    ),
    number
    (
        "n",
        "number",

        "Parameter n that is specific for each polyhedra class.",

        false,
        0,
        "integer number",
        cmd
    ),
    random
    (
        "",
        "random",

        "Switch generation from concrete types of polyhedrons to random vertices/inequalities generation inside polyhedron or bound it respectively. "
        "An argument is an positive integer number that denotes the number of elements in the random set (number of vertices, inequalities etc. -- "
        "all depends on class of polyhedra).",

        false,
        0,
        "positive integer number",
        cmd
    ),
    repulsion
    (
        "",
        "repulsion",

        "Applied to randomly generated points. Applay several iterations of repulsion process to points making their distribution more regular in the space. "
        "Number provided is the number of iteration of repulsion process. Value 0 means repultion turned off. Currently applicable to sphere type only.",

        false,
        0,
        "positive integer number",
        cmd
    )
    /*,
    vertices
    (
        "",
        "vertices",

        "Switch output representation to vertices. For some kind of polyhedra it can be very expensive.",

        cmd,
        false
    ),
    facets
    (
        "",
        "facets",

        "Switch output representation to facets. For some kind of polyhedra it can be very expensive.",

        cmd,
        false
    )*/
    {
    }


}
}
}
