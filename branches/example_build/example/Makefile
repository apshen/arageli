##############################################################################
#
#    Makefile
#
#    This makefile builds examples of arageli using
#
#    This file is a part of the Arageli library.
#
#    Copyright (C) 2006--2007 Aleksey Bader
#    University of Nizhni Novgorod, Russia
#
#    The Arageli Library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License version 2
#    as published by the Free Software Foundation.
#
#    The Arageli Library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
#
#    We are also open for dual licensing for the whole library or
#    for its particular part. If you are interested to get the library
#    in this way, i.e. not under the GNU General Public License,
#    please contact Arageli Support Service support.arageli@gmail.com.
#
##############################################################################

BIN=../bin
INCLUDE_DIR=../src
LIBRARY_DIR=../lib

all: \
     $(BIN)/big_int_basic\
     $(BIN)/big_int_rational_creation\
     $(BIN)/big_int_rational_functions\
     $(BIN)/big_int_rational_operations\
     $(BIN)/adv_simplex_method\

$(BIN)/adv_simplex_method: adv_simplex_method.cpp adv_simplex_method.task.txt
	$(CXX) -I"$(INCLUDE_DIR)" adv_simplex_method.cpp -o$@ -L"$(LIBRARY_DIR)" -larageli

$(BIN)/big_int_basic: big_int_basic.cpp
	$(CXX) -I"$(INCLUDE_DIR)" big_int_basic.cpp -o$@ -L"$(LIBRARY_DIR)" -larageli

$(BIN)/big_int_rational_creation: big_int_rational_creation.cpp big_int_rational_creation.cpp.in.txt
	$(CXX) -I"$(INCLUDE_DIR)" big_int_rational_creation.cpp -o$@ -L"$(LIBRARY_DIR)" -larageli

$(BIN)/big_int_rational_functions: big_int_rational_functions.cpp
	$(CXX) -I"$(INCLUDE_DIR)" big_int_rational_functions.cpp -o$@ -L"$(LIBRARY_DIR)" -larageli

$(BIN)/big_int_rational_operations: big_int_rational_operations.cpp
	$(CXX) -I"$(INCLUDE_DIR)" big_int_rational_operations.cpp -o$@ -L"$(LIBRARY_DIR)" -larageli
