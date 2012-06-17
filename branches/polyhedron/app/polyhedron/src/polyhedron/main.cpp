/*****************************************************************************

    main.cpp

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

#include "ball.hpp"
#include "CC2d.hpp"
#include "cube.hpp"
#include "cyclic.hpp"
#include "dodecahedron.hpp"
#include "DPP2d.hpp"
#include "dwarfed_cube.hpp"
#include "hecatonicosachoron.hpp"
#include "hexacosichoron.hpp"
#include "icosahedron.hpp"
#include "icositetrachoron.hpp"
#include "octahedron.hpp"
#include "simplex.hpp"
#include "sphere.hpp"
#include "TT2d.hpp"
#include "transport.hpp"

using namespace Arageli::app::polyhedron;

int main (int argc, char** argv)
{
    try
    {
        // Fill the table with acceptable types of polyhedrons.
        typedef std::map<std::string, Processor*> Types;
        std::map<std::string, Processor*> types;
        types["simplex"] = new Simplex;
        types["tetrahedron"] = types["simplex"];
        types["cube01"] = new Cube01;
        types["cyclic"] = new Cyclic;
        types["dwarfedcube"] = new DwarfedCube;
        types["dwarfcube"] = types["dwarfedcube"];
        types["CC2d"] = new CC2d;
        types["DPP2d"] = new DPP2d;
        types["TT2d"] = new TT2d;
        types["cube"] = new Cube;
        types["octahedron"] = new Octahedron;
        types["dodecahedron"] = new Dodecahedron;
        types["icosahedron"] = new Icosahedron;
        types["icositetrachoron"] = new Icositetrachoron;
        types["hecatonicosachoron"] = new Hecatonicosachoron;
        types["hexacosichoron"] = new Hexacosichoron;
        types["ball"] = new Ball;
        types["sphere"] = new Sphere;
        types["transport"] = new Transport;

        TCLAP::CmdLine cmd
            (
            "You are running Polyhedron 0.1 application.\n"
            "Polyhedron is a command line interface to "
            "a set of polyhedron generators.\n"
            "Copyright (C) Sergey S. Lyalin, 2009-2010\n"
            "Copyright (C) Anastasya A. Ryzhova, 2010\n",

            ' ',

            "0.1.1"
            );

        CmdArgs cmdargs(cmd);

        TCLAP::ValueArg<std::string> output
            (
            "o",
            "output",

            "Name of output file.",

            false,
            "polyhedron.out",
            "acceptable file name",
            cmd
            );

        cmd.parse(argc, argv);

        Types::iterator type_iter = types.find(cmdargs.type.getValue());
        if(type_iter != types.end())
        {
            std::ofstream file(output.getValue().c_str());
            type_iter->second->generate(file, cmdargs);
        }
        else
        {
            std::cerr << "[ERROR] Unknown type `" << cmdargs.type.getValue() << "'.";
        }

        // WARNING! MEMORY LEAKS FOR types.
    }
    catch(const char* msg)
    {
        std::cerr << "[ ERROR ] " << msg << '\n';
        return 2;
    }
    catch(...)
    {
        std::cerr
            << "\nINTERNAL ERROR\n"
            << "Unknown exception has been thrown.\n";
        return 1;
    }
}
