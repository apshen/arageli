/*****************************************************************************

    main.cpp

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

#include "compare_polyhedra.hpp"

using namespace Arageli::app::phdiff;

int main (int argc, char** argv)
{
    try
    {
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

        TCLAP::ValueArg<std::string> input1
        (
            "",
            "input1",

            "Name of the first input file.",

            false,
            "polyhedron1.in",
            "acceptable file name",
            cmd
        );

		TCLAP::ValueArg<std::string> input2
        (
            "",
            "input2",

            "Name of the second input file.",

            false,
            "polyhedron2.in",
            "acceptable file name",
            cmd
        );

        cmd.parse(argc, argv);

        if( cmdargs.type.getValue() == "int" || cmdargs.type.getValue() == "bigint" || cmdargs.type.getValue() == "double")
        {
            std::ifstream file1(input1.getValue().c_str());
			std::ifstream file2(input2.getValue().c_str());
            if ( compare_polyhedra(file1, file2, cmdargs) )
			{
				std::cerr << "Polyhedra are identical.\n"; 
			}
			else
			{
				std::cerr << "Polyhedra are not identical.\n";
			}

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
    catch(const cannot_find_split_dim&)
    {
        std::cerr << "[ ERROR ] Cannot find split; the first input matrix is ill-formed.\n";
    }
    catch(...)
    {
        std::cerr
            << "\nINTERNAL ERROR\n"
            << "Unknown exception has been thrown.\n";
        return 1;
    }
}
