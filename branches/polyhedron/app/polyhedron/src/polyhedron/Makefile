##############################################################################
#
#    Makefile
#
#    Makefile for GNU Make Utility to build the library and
#    other programs from the Arageli Package and run some test and
#    examples.
#
#    This file is a part of the Arageli library.
#
#    Copyright (C) 2010 Anastasya A. Ryzhova
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

ARAGELI_DIR = ../../../..
SOURCE_DIR = .

CXXFLAGS += -MD -g -DARAGELI_SOURCE -O2 -DARAGELI_DEBUG_LEVEL=3 -DARAGELI_ASSERT_THROW_EXCEPTION
LINKFLAGS = -L"$(ARAGELI_DIR)/lib/" -ltestarageli

all: $(ARAGELI_DIR)/bin/polyhedron

$(ARAGELI_DIR)/bin/polyhedron: $(ARAGELI_DIR)/lib/libtestarageli.a \
	stdafx.hpp.gch \
	$(patsubst %.cpp,%.o, $(wildcard $(addsuffix /*.cpp,$(SOURCE_DIR))))
	make -C $(ARAGELI_DIR)/src/arageli ../../lib/libtestarageli.a
	$(CXX) -o $(ARAGELI_DIR)/bin/polyhedron $(patsubst %.cpp,%.o, $(wildcard $(addsuffix /*.cpp,$(SOURCE_DIR)))) $(LINKFLAGS)

stdafx.hpp.gch: stdafx.hpp
	$(CC) -x c++-header $(CXXFLAGS) -o $@ $< -I"$(ARAGELI_DIR)/src/" -I"$(ARAGELI_DIR)/tool/tclap/include/"

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ -I"$(ARAGELI_DIR)/src/" -I"$(ARAGELI_DIR)/tool/tclap/include/"

$(ARAGELI_DIR)/lib/libtestarageli.a:
	$(MAKE) -C $(ARAGELI_DIR)/src/arageli ../../lib/libtestarageli.a

TAGS:
	ctags $(wildcard $(addsuffix /*.cpp,$(SOURCE_DIR))) $(wildcard $(addsuffix /*.hpp,$(SOURCE_DIR)))

clean:
	rm -f $(ARAGELI_DIR)/bin/polyhedron
	rm -f *.o *.d
	
cleanall: clean
	rm -f stdafx.hpp.gch
	rm -f tags

include $(wildcard *.d)
