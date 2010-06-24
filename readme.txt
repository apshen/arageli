********************************************************************

                The Arageli Library Readme File

********************************************************************

Copyright (C) Nikolai Yu. Zolotykh, 1999--2009
Copyright (C) Sergey S. Lyalin, 2005--2009
University of Nizhni Novgorod, Russia.

Arageli is a C++ library and a package of programs for
computations in arithmetic, algebra, geometry, linear and integer
linear programming. Arageli is a library for dealing with precise,
i.e. symbolic or algebraic, computations. It contains a definition
to model basic algebraic structures such as integer numbers with
arbitrary precision, rational numbers, vectors, matrices,
polynomials etc. Arageli is written in C++ and use power and
expressiveness of the language.

The home page: http://www.arageli.org

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

********************************************************************


1. Software Requirements

To build the library and create documentation you need several
applications installed. Without some of them you will not
able to create some parts of the documentation or the library
itself.

1.1. To compile the library

You need one of the C++ 1998 standard compilers. We have checked
compilation with:
    - GCC 4.1.3 20070929 (prerelease) (Ubuntu 4.1.2-16ubuntu2),
    - Microsoft (R) 32-bit C/C++ Optimizing Compiler Version 15
      for 80x86 (included in Visual Studio 2008),
    - Microsoft (R) C/C++ Optimizing Compiler Version 15 for x64
      (included in Visual Studio 2008).
Please let us know if you faced troubles with those compilers or
you have success or faulure to build with other C++ compilers/platforms
to include or exclude them to/from the list.

1.2 To create all guides (not auto generated)

You need LaTeX installation. We use MiKTeX and teTex.
Note that some documentation files are written in Russian.

1.3 To create guides with the use of lgrind

You need lgrind installation. Without this you cannot create
some guide files with the lgrind directives.

1.4 To create a reference to sources (auto generated)

You need doxygen installation. Without this you cannot extract
documented items from sources.

1.5 To create Arageli User's Guide

You need to be complied with requirements from 1.2, 1.3, 1.4;
Microsoft Windows environment. Without Microsoft Windows Shell
you need to compile a lot of files manually.

To create all parts of the distributive you should have all
programs referenced at 1.1--1.4. If you plan to contribute into
the library, we recommend you to have all these programs and
additionally SubVersion client installation to be able to commit
changes to the repository.

If you need some part of the library compiled (both library or
documentation), you can download it at the home page of Arageli
http://www.arageli.org.


2. Building

Almost all build scripts are located in build directory of
the distributive but some of them are spread along the directories.
Choose proper platform directory for your system and build all
or a particular part of the distributive.

If you use Microsoft Visual Studio 2005, you can try to build
solution that located at build/msvs_2005/arageli_all.sln. Please
be informed that support of Visual Studio 2005 is deprecated.

If you use Microsoft Visual Studio 2008, you can build
solution that located at build/msvs_2008/arageli_all.sln.

If you are in some Linux system, you can build the library,
test system and test set by the make utility. To make only
Arageli library go to library root directory and type:

    make

After that, libarageli.a file appears in the lib directory
of the package. To make test set, type:

    make check

After that, test executive file appears in the bin directory
of the package. Note, while executing the last command
arageli and ts systems have being built. To make tests and
run them, type:

    make runtests

After that, test.log file appears in the status directory
of the package.

You can use CXXFLAGS, ARAGELICXXFLAGS and TESTSCXXFLAGS
variables to set specific command line parameters to tune
compilation of different parts as usual.


3. Using

To use the library in your application you have to include
one or several necessary header files from ./src/arageli
directory. If you do not know which files you need to use
a particulary feature of Arageli, you can just inlude
./src/arageli/arageli.hpp that includes all header files of
the library.

If you use Arageli in an executable application (not a static
library) or in a dynamic library, you have to link your
target binary with Arageli static library placed in ./lib
directory. There are several Arageli static libraries;
see parameters of the build script to choose right version of
static library file to link with your application/library.


4. Additional Documentation and Feedback

Some ready-to-use documentation files already present in doc
directory of the distributive. Go to http://www.arageli.org
for updates and news. Mail questions, comments and suggestions on
support.arageli@gmail.com.
